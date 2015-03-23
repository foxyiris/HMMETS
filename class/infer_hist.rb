class InferHist
  include Math

  attr_reader :gain_branch
  attr_reader :signal_gain_branch
  attr_reader :gene_lost_branch
  attr_reader :signal_lost_branch

  MAX     = 1e50
  FP_rate = 0.017
  FN_rate = 0.2

  def initialize(pt = phylo_tree, pp = phylo_profile, input_genes)
    @pt          = pt
    @pp          = pp
    @N           = 2*@pt.n_s-1
    @input_genes = input_genes   # this is an array which stores the names of the target genes
    @is_clade    = Array.new(@N).map{ Array.new(@N) }
    @null_m      = [
                    [log(1), log(0)],
                    [log(0.03), log(0.97)]
                   ]
    # below var contains sampled states at each node for each gene
    @hidden_samp = Array.new(input_genes.size).map{ Array.new(@N, nil) }

    # parameters for prior
    @a   = 0.0045 # beta
    @b   = 0.1455 # beta
    @al1 = 0.001  # dirichlet
    @al2 = 0.002  # dirichlet
    @al3 = 0.003  # dirichlet

    # in the future this has to be k-D array? to store MULTI gainings.
    @gain_branch        = Hash.new()
    @signal_gain_branch = Hash.new()
    # lost branch can be many, so they are contained as an array
    @gene_lost_branch   = Hash.new {|hash,key| hash[key] = Array.new()}
    @signal_lost_branch = Hash.new {|hash,key| hash[key] = Array.new()}

  end

  def pre_process
    # estimate gain of one property
    @input_genes.each do |gene|
      calc_gene_gain_org(gene)
      calc_signal_gain_org(gene)
      #joint_ML(gene)
      #maximum_parsimony(gene)
    end
  end

  def mcmc(total, burn_in)
    (0..total-1).each do |itr|
      STDERR.puts "Running..#{itr}"
      @input_genes.each_with_index do |gene, i|
        sampling_state(i, gene)
      end
    end
  end

  def marginal_ML(g = gene, gn = g_gain_node, sn = s_gain_node)

    out_of_gg_clade = 0
    out_of_sg_clade = 0

    # memo for me-> Ruby assings same reference with such declarations:
    #   Array.new(N, Array.new(K,0)). In this case, N arrays have the same reference.
    #   Confirm not to declare with above style in ruby when using multi-D array.
    log_prob_from_leaves = Array.new(@N).map{ Array.new(2).map{ Array.new(3, 0) }}
    log_prob_near_leaves = Array.new(@N).map{ Array.new(3,0) }
    log_prob             = Array.new(@N).map{ Array.new(3,0) }

    #sg    = 0.0072 # signal gain, global
    #sg    = 0.0830 # signal gain, mts annotated
    sl    = 0.07554 # signal loss, global
    #sl    = 0.05692 # signal loss, mts annotated
    gl    = 0.02    # gene loss
    eps   = 0.01
    eps_m = 0.08

    prob_pred_error = [
                       [log(1-eps), log(eps/2),         log(eps/2)],
                       [log(eps/2), log(1-eps_m-eps/2), log(eps_m)],
                       [log(eps/2), log(eps_m),         log(1-eps_m-eps/2)]
                      ]

    t_m = [
           [log(1),  log(0),    log(0)],
           [log(gl), log(1-gl), log(0)],
           [log(gl), log(sl),   log(1-sl-gl)]
          ]

    # quantitate fluctuation of observed states by pre-calculated performance.
    (0..@pt.n_s-1).each{ |num|

      node  = @pt.sorted_nodes[num]
      state = @pp.profiles[@pp.symbol2index[g]][num]

      if gn != node && !@pt.tree.descendents(gn, root=@pt.root).include?(node)
        out_of_gg_clade += prob_pred_error[0][state]
        next
      elsif gn == node || @pt.tree.descendents(gn, root=@pt.root).include?(node)
        if sn != node && !@pt.tree.descendents(sn, root=@pt.root).include?(node)
          out_of_sg_clade += prob_pred_error[1][state]
          next
        end
      end

      log_prob_near_leaves[num][0] = prob_pred_error[0][state]
      log_prob_near_leaves[num][1] = prob_pred_error[1][state]
      log_prob_near_leaves[num][2] = prob_pred_error[2][state]
    }

    n_index = @pt.sorted_nodes.index(sn)
    @pt.sorted_nodes.each_with_index { |node, num|
      if num >= n_index
        break
      end

      if gn != node && !@pt.tree.descendents(gn, root=@pt.root).include?(node)
        next
      elsif gn == node || @pt.tree.descendents(gn, root=@pt.root).include?(node)
        if sn != node && !@pt.tree.descendents(sn, root=@pt.root).include?(node)
          next
        end
      end

      log_prob[num][0] = log_prob_near_leaves[num][0] + log_prob_from_leaves[num][0][0] + log_prob_from_leaves[num][1][0]
      log_prob[num][1] = log_prob_near_leaves[num][1] + log_prob_from_leaves[num][0][1] + log_prob_from_leaves[num][1][1]
      log_prob[num][2] = log_prob_near_leaves[num][2] + log_prob_from_leaves[num][0][2] + log_prob_from_leaves[num][1][2]

      anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]
      chd_index = @pt.which_child[node]

      log_prob_from_leaves[anc_index][chd_index][0] = Utils.log_sum_exp3(
                                                                        t_m[0][0]+log_prob[num][0],
                                                                        t_m[0][1]+log_prob[num][1],
                                                                        t_m[0][2]+log_prob[num][2])

      log_prob_from_leaves[anc_index][chd_index][1] = Utils.log_sum_exp3(
                                                                        t_m[1][0]+log_prob[num][0],
                                                                        t_m[1][1]+log_prob[num][1],
                                                                        t_m[1][2]+log_prob[num][2])

      log_prob_from_leaves[anc_index][chd_index][2] = Utils.log_sum_exp3(
                                                                        t_m[2][0]+log_prob[num][0],
                                                                        t_m[2][1]+log_prob[num][1],
                                                                        t_m[2][2]+log_prob[num][2])
    }

    log_prob[n_index][0] = log_prob_from_leaves[n_index][0][0] + log_prob_from_leaves[n_index][1][0]
    log_prob[n_index][1] = log_prob_from_leaves[n_index][0][1] + log_prob_from_leaves[n_index][1][1]
    log_prob[n_index][2] = log_prob_from_leaves[n_index][0][2] + log_prob_from_leaves[n_index][1][2]

    return log_prob[n_index][2] + out_of_gg_clade + out_of_sg_clade
  end # end of method


  # 2015/3/12 Pupko 2000
  # In the field of ancestral reconstruction, below maximum parsimony method is nowdays old-fashioned.
  # I try to reconstruct ancestor state of signal with ML and transition parameter determined by below function.
  # But problem is that we still ignore branch specific transition or length, so bayesian way should be an appropriate.
  # In this case, reconstruction has to be done with many matrices, so this method cannot be applied to the model.
  # In such a case, Pupko 2002 method should be better in terms of assumption and its speed with their approximation.

  # 2015/3/13
  # I realized joint ML method cannot work with the transition matix with initial assumption.
  # Since t00 is too huge, namely 1, root is likely to be 0 and descendents are also likely to be 0 due to the prob.
  # As a result, even though I gave a transition matrix, this violates observation.
  # maybe this is why CLIME assumes one gaining and it doesn't condider L(0) to pick up gain branch. (it is easy L(0) > L(1), and max will be 0) 
  def joint_ML(g=gene)

    # considering gene loss model but assume no regaining of a gene
    
    #sg    = 0.0072 # signal gain, global
    sg    = 0.0830 # signal gain, mts annotated
    #sl    = 0.07554 # signal loss, global
    sl    = 0.05692 # signal loss, mts annotated
    gl    = 0.02    # gene loss
    eps   = 0.01
    eps_m = 0.08

    prob_pred_error = [
                       [log(1-eps), log(eps/2),         log(eps/2)],
                       [log(eps/2), log(1-eps_m-eps/2), log(eps_m)],
                       [log(eps/2), log(eps_m),         log(1-eps_m-eps/2)]
                      ]

    t_m = [
           [log(1),  log(0),       log(0)],
           [log(gl), log(1-sg-gl), log(sg)],
           [log(gl), log(sl),      log(1-sl-gl)]
          ]

    # this is pre-calculated.
    gain_node_num = @pt.node2num[@gain_branch[g]]

    # counter for transition event
    counter = Array.new(3).map{ Array.new(3,0) }

    #each node has cost vector where vector[i] is minimum cost when nodes' parent is state i.
    #S_i(x), namely contain minimum cost when node i's parent is x
    likel_vec             = Array.new(@N).map{ Array.new(3,0) }
    #C_i(x), namely contain i's state giving minimum cost when parent is x
    #element can be multiple indices.
    most_likel_state_vec  = Array.new(@N).map{ Array.new(3) }
    
    likel_vec_from_leaves = Array.new(@N).map{ Array.new(2).map{ Array.new(3,0) } }

    # calc cost vector for leaves
    @pt.sorted_nodes.each_with_index {|node, num|

      if num > @pt.n_s-1
        break
      end

      # we don't believe even leaf state. complete prob model.
      # but in this case, this is the same as observation is correct premice after all.
      # maybe we have to take forward-summation-backward-sampling like method.
      state_at_leaf = @pp.profiles[@pp.symbol2index[g]][num]
      (0..2).each do |state| #leaf parents is state
        max       = log(0)
        max_state = -1
        (0..2).each do |y|
          prob = prob_pred_error[state][y]
          #puts "debug: #{num} #{state} #{y} #{prob} #{max} #{max_state}"
          if prob > max
            max       = prob
            max_state = y
          end
        end

        likel_vec[num][state] = max
        most_likel_state_vec[num][state] = max_state
      end

      #likel_vec[num][0] = prob_pred_error[0][state_at_leaf]
      #likel_vec[num][1] = prob_pred_error[1][state_at_leaf]
      #likel_vec[num][2] = prob_pred_error[2][state_at_leaf]

      anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]
      chd_index = @pt.which_child[node]

      # state at leaf is always observation after all even we take prob error into account.
      # max_i P(H_{g,leaf}=i}|Xg,s) = observation state
      #most_likel_state_vec[num][0] = [state_at_leaf]
      #most_likel_state_vec[num][1] = [state_at_leaf]
      #most_likel_state_vec[num][2] = [state_at_leaf]

      likel_vec_from_leaves[anc_index][chd_index][0] = likel_vec[num][0]
      likel_vec_from_leaves[anc_index][chd_index][1] = likel_vec[num][1]
      likel_vec_from_leaves[anc_index][chd_index][2] = likel_vec[num][2]

    }

    # then calc from the bottom for internal node other than root
    @pt.sorted_nodes.each_with_index { |node, num|

      if num <= @pt.n_s-1
        next
      elsif num  == gain_node_num
        break # this calculation is conditioned by lambda
      elsif node == @pt.root
        break
      end

      (0..2).each do |state|

        # assume there is only one optimum state.
        max       = log(0) # since cost is increasing, MAX is not enough big to choose min. this is an awkward way tho. 
        max_state = -1
        # S_i(x) = min_y [c(x,y) + S_left_child(y) + S_right_child(y)]
        (0..2).each do |y|
          prob = t_m[state][y] + likel_vec_from_leaves[num][0][y] + likel_vec_from_leaves[num][1][y]
          #puts "debug: #{num} #{state} #{y} #{prob} #{max} #{max_state}"
          if prob > max
            max = prob
            max_state = y
          end
        end

        likel_vec[num][state] = max.to_f
        if max_state == -1
          STDERR.puts "Error: max state was not updated for #{state} at node #{num}"
          p likel_vec[num]
        else
          most_likel_state_vec[num][state] = max_state
        end


      end


      anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]
      chd_index = @pt.which_child[node]

      likel_vec_from_leaves[anc_index][chd_index][0] = likel_vec[num][0]
      likel_vec_from_leaves[anc_index][chd_index][1] = likel_vec[num][1]
      likel_vec_from_leaves[anc_index][chd_index][2] = likel_vec[num][2]

    }

    root_num   = gain_node_num # again, this ML is conditioned by lambda. root_num can be either root or descendents of root.
    root_state = -1
    max        = log(0)

    # since root_num is conditioned with state 1 or 2. we can ignore 0.
    (1..2).each do |state|
      #puts "state: #{state} #{cost_vec_from_leaves[root_num][0][state]} #{cost_vec_from_leaves[root_num][1][state]}"
      likel = likel_vec_from_leaves[root_num][0][state] + likel_vec_from_leaves[root_num][1][state]
      if likel > max
        max        = likel
        root_state = state.to_i
      end

      likel_vec[root_num][state] = likel
    end

    puts "root #{root_state} #{likel_vec[root_num].join(',')}"

    @pt.tree.children(@gain_branch[g], root=@pt.root).each do |child|
      _trace_back_state(child, most_likel_state_vec, root_state, counter)
    end

    #p most_likel_state_vec

  end

  def _trace_back_state(n=node, s_vec, state, counter)

    # state is scalar
    array = @pt.tree.children(n, root=@pt.root)
    if array.size != 0
      array.each do |child|
        #counter[candidate][s_vec[@pt.node2num[child]][candidate][0]] += weight
        #puts "#{@pt.node2num[n]}:#{state} -> #{@pt.node2num[child]}:#{s_vec[@pt.node2num[child]][state]}"
        _trace_back_state(child, s_vec, s_vec[@pt.node2num[child]][state], counter)
      end
    end

  end


  # 2015/3/10
  # I faced a problem to estimate parameters of null probs, namely P01 and P10 (1=MTS, 0=Non-MTS)
  # Above null model apperantly has no objective reasoning (P10=0.03 and P01=0).
  # To avoid theoretical problems, I have to start some principle based method like maximum parsimony to get some flavor.
  # (No experimental or statistical report to argue P01 and P10, so I have to start from the beginning.)
  # To calculate fast, DP is applied with a global cost matrix Transition (01 and 10) has 1 otherwise 0.
  # If an ancestral node has equal cost for both states, what should I do? At present, count up both with half weight, namely 0.5.

  # 2015/3/11
  # With a condition without any gene loss in 20 yeast speacies, P01=0.00759 and P10=0.0711
  # Since MitoFates's fn rate is a bit high, so P10 might be overestimated, but anyway this is a start line.
  def maximum_parsimony(g = gene)
    # without gene loss model
    #cost_m = [
    #          [0,1],
    #          [1,0]]

    # considering gene loss model but assume no regaining of a gene
    cost_m = [
              [0,MAX,MAX],
              [1,0,1],
              [1,1,0]]

    # considering gene loss model without any assumption
    #cost_m = [
    #          [0,1,1],
    #          [1,0,1],
    #          [1,1,0]]

    # counter for transition event
    counter = Array.new(3).map{ Array.new(3,0) }

    #each node has cost vector where vector[i] is minimum cost when nodes' parent is state i.
    #S_i(x), namely contain minimum cost when node i's parent is x
    cost_vec             = Array.new(@N).map{ Array.new(3,0) }
    #C_i(x), namely contain i's state giving minimum cost when parent is x
    #element can be multiple indices.
    min_cost_state_vec   = Array.new(@N).map{ Array.new(3) }
    
    cost_vec_from_leaves = Array.new(@N).map{ Array.new(2).map{ Array.new(3,0) } }

    # calc cost vector for leaves
    @pt.sorted_nodes.each_with_index {|node, num|

      if num > @pt.n_s-1
        break
      end

      state_at_leaf = @pp.profiles[@pp.symbol2index[g]][num]
      cost_vec[num][0] = cost_m[0][state_at_leaf]
      cost_vec[num][1] = cost_m[1][state_at_leaf]
      cost_vec[num][2] = cost_m[2][state_at_leaf]

      anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]
      chd_index = @pt.which_child[node]

      # state at leaf is always observation.
      min_cost_state_vec[num][0] = [state_at_leaf]
      min_cost_state_vec[num][1] = [state_at_leaf]
      min_cost_state_vec[num][2] = [state_at_leaf]

      cost_vec_from_leaves[anc_index][chd_index][0] = cost_vec[num][0]
      cost_vec_from_leaves[anc_index][chd_index][1] = cost_vec[num][1]
      cost_vec_from_leaves[anc_index][chd_index][2] = cost_vec[num][2]

    }

    # then calc from the bottom for internal node other than root
    @pt.sorted_nodes.each_with_index { |node, num|

      if num <= @pt.n_s-1
        next
      elsif node == @pt.root
        break
      end

      (0..2).each do |state|

        min       = MAX*MAX # since cost is increasing, MAX is not enough big to choose min. this is an awkward way tho. 
        min_state = Array.new()
        min_state[0] = -1
        # S_i(x) = min_y [c(x,y) + S_left_child(y) + S_right_child(y)]
        (0..2).each do |y|
          cost = cost_m[state][y] + cost_vec_from_leaves[num][0][y] + cost_vec_from_leaves[num][1][y]
          #puts "debug: #{num} #{state} #{y} #{cost}"
          if cost < min
            min = cost
            min_state[0] = y
          elsif cost == min
            min_state.push(y)
          end
        end

        cost_vec[num][state] = min.to_i
        if min_state[0] == -1
          STDERR.puts "Error: min state was not updated at node #{num}"
          p cost_vec[num]
        elsif min_state.size != 1
          if min_state.uniq.size == 1
            STDERR.puts "Error: multi min state is homogeneous with #{min_state[0]} at node #{num}"
          end
          min_cost_state_vec[num][state] = min_state
        elsif min_state.size == 1
          min_cost_state_vec[num][state] = min_state
        end


      end


      anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]
      chd_index = @pt.which_child[node]

      cost_vec_from_leaves[anc_index][chd_index][0] = cost_vec[num][0]
      cost_vec_from_leaves[anc_index][chd_index][1] = cost_vec[num][1]
      cost_vec_from_leaves[anc_index][chd_index][2] = cost_vec[num][2]

    }

    root_num = @pt.node2num[@pt.root]

    (0..2).each do |state|
      #puts "state: #{state} #{cost_vec_from_leaves[root_num][0][state]} #{cost_vec_from_leaves[root_num][1][state]}"
      cost = cost_vec_from_leaves[root_num][0][state] + cost_vec_from_leaves[root_num][1][state]
      cost_vec[root_num][state] = cost.to_i
    end
    
    min_state = cost_vec[root_num].each.with_index.map{ |a, i| (a == cost_vec[@pt.node2num[@pt.root]].min) ? i : nil }.compact.uniq
    @pt.tree.children(@pt.root, root=@pt.root).each do |child|
      _get_cost_and_state(child, min_cost_state_vec, min_state, 1.0, counter)
    end

    puts counter.join(", ")
    
  end

  def _get_cost_and_state(n=node, s_vec, state, weight, counter)

    state.each do |candidate|
      array = @pt.tree.children(n, root=@pt.root)
      if array.size != 0
        array.each do |child|
          #puts "#{@pt.node2num[n]}:#{candidate} -> #{@pt.node2num[child]}:#{s_vec[@pt.node2num[child]][candidate][0]}, #{weight}"
          counter[candidate][s_vec[@pt.node2num[child]][candidate][0]] += weight
          _get_cost_and_state(child, s_vec, s_vec[@pt.node2num[child]][candidate], weight/state.size, counter)
        end
      end
    end

  end

  def calc_gene_gain_org(g = gene)
    max = log(0)
    gain_node = nil

    # tentatively, assume one gaining.
    # loop with all nodes.
    @pt.sorted_nodes.each do |node|
      prob = calc_log_Xg(g, node)
      #STDERR.puts "DEBUG: #{g} #{prob}"
      if prob > max
        max = prob
        gain_node = node
      end
    end
    
    puts "MAX log likelihood: #{max} for #{g}"
    @gain_branch[g] = gain_node

  end

  def calc_signal_gain_org(g = gene)
    max = log(0)
    signal_gain_node = nil
    gene_gain_node   = @gain_branch[g]
    # tentatively, assume one gaining.
    # loop with all nodes.
    @pt.sorted_nodes.each do |node|
      if gene_gain_node == node || @pt.tree.descendents(gene_gain_node, root=@pt.root).include?(node)
        prob = marginal_ML(g, @gain_branch[g], node)
        #STDERR.puts "DEBUG: #{g} #{prob}"
        if prob > max
          max = prob
          signal_gain_node = node
        end
      end
    end
    
    puts "MAX log likelihood given gene gain node: #{max} for #{g}"
    @signal_gain_branch[g] = signal_gain_node

  end

  def sampling_state(i = input_index, g = gene)
    gg = @gain_branch[g]
    sg = @signal_gain_branch[g]

    if !gg || !sg
      STDERR.puts "Error: gene gain or/and signal gain have not been estimated."
      return
    end

    # memo for me-> Ruby assings same reference with such declarations:
    #   Array.new(N, Array.new(K,0)). In this case, N arrays have the same reference.
    #   Confirm not to declare with above style in ruby when using multi-D array.
    log_prob_from_leaves     = Array.new(@N).map{ Array.new(2).map{ Array.new(3, 0) }}
    log_prob_near_leaves     = Array.new(@N).map{ Array.new(3,0) }
    log_prob                 = Array.new(@N).map{ Array.new(3,0) }
    log_backward             = Array.new(@N).map{ Array.new(3,0) }
    log_backward_from_parent = Array.new(@N).map{ Array.new(3,0) }

    # below var contains branch specific transition params for each input gene
    t_mat = Array.new(@N).map{
              Array.new(3).map{
                Array.new(3, log(0))
              }
            }


    #sg    = 0.0072 # signal gain, global
    #sg    = 0.0830 # signal gain, mts annotated
    sl    = 0.07554 # signal loss, global
    #sl    = 0.05692 # signal loss, mts annotated
    gl    = 0.02    # gene loss
    eps   = 0.01
    eps_m = 0.08

    prob_pred_error = [
                       [log(1-eps), log(eps/2),         log(eps/2)],
                       [log(eps/2), log(1-eps_m-eps/2), log(eps_m)],
                       [log(eps/2), log(eps_m),         log(1-eps_m-eps/2)]
                      ]

    t_m = [
           [log(1),  log(0),    log(0)],
           [log(gl), log(1-gl), log(0)],
           [log(gl), log(sl),   log(1-sl-gl)]
          ]

    if @hidden_samp[i][0] == nil
    #if true
      # in case of first round
      (0..@N-1).each do |i|
        t_mat[i] = t_m
      end
    else
      # estimate branch specific parameter
      @pt.sorted_nodes.each_with_index { |node, num|
        if num >= @N-1
          break
        end

        if gg != node && !@pt.tree.descendents(gg, root=@pt.root).include?(node)
          next
        elsif gg == node || @pt.tree.descendents(gg, root=@pt.root).include?(node)
          if sg != node && !@pt.tree.descendents(sg, root=@pt.root).include?(node)
            # non-mts region
            count = Array.new(2).map{ Array.new(2,0) }
            anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]

            count[@hidden_samp[i][anc_index]][@hidden_samp[i][num]] += 1

            t_mat[num][0][0] = 0
            t_mat[num][1][0] = log((@a+count[1][0])/(@a+@b+count[1][0]+count[1][1])) #beta dist
            t_mat[num][1][1] = log(1 - t_mat[num][1][0])
            t_mat[num][1][2] = log(0)
            t_mat[num][2][0] = log(0)
            t_mat[num][2][1] = log(0)
            t_mat[num][2][2] = log(0)

          else
            # mts region
            count = Array.new(3).map{ Array.new(3,0) }
            anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]

            count[@hidden_samp[i][anc_index]][@hidden_samp[i][num]] += 1

            t_mat[num][0][0] = 0
            t_mat[num][1][0] = log((@a+count[1][0])/(@a+@b+count[1][0]+count[1][1])) #beta dist
            t_mat[num][1][1] = log(1 - exp(t_mat[num][1][0]))
            t_mat[num][1][2] = log(0)
            t_mat[num][2][0] = log((@al1+count[2][0])/(@al1+@al2+@al3+count[2][0]+count[2][1]+count[2][2])) #gamma dist
            t_mat[num][2][1] = log((@al2+count[2][1])/(@al1+@al2+@al3+count[2][0]+count[2][1]+count[2][2]))
            t_mat[num][2][2] = log(1 - exp(t_mat[num][2][0]) - exp(t_mat[num][2][1]))

          end
        end
      } # end of each

    end # end of if

    # quantitate fluctuation of observed states by pre-calculated performance.
    (0..@pt.n_s-1).each{ |num|

      node  = @pt.sorted_nodes[num]
      state = @pp.profiles[@pp.symbol2index[g]][num]

      log_prob_near_leaves[num][0] = prob_pred_error[0][state]
      log_prob_near_leaves[num][1] = prob_pred_error[1][state]
      log_prob_near_leaves[num][2] = prob_pred_error[2][state]
    }

    n_index = @pt.sorted_nodes.index(gg)
    s_index = @pt.sorted_nodes.index(sg)

    @pt.sorted_nodes.each_with_index { |node, num|
      if num >= n_index
        break
      end

      if gg != node && !@pt.tree.descendents(gg, root=@pt.root).include?(node)
        next
      elsif gg == node || @pt.tree.descendents(gg, root=@pt.root).include?(node)
        if sg != node && !@pt.tree.descendents(sg, root=@pt.root).include?(node)
          next
        end
      end

      log_prob[num][0] = log_prob_near_leaves[num][0] + log_prob_from_leaves[num][0][0] + log_prob_from_leaves[num][1][0]
      log_prob[num][1] = log_prob_near_leaves[num][1] + log_prob_from_leaves[num][0][1] + log_prob_from_leaves[num][1][1]
      log_prob[num][2] = log_prob_near_leaves[num][2] + log_prob_from_leaves[num][0][2] + log_prob_from_leaves[num][1][2]

      anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]
      chd_index = @pt.which_child[node]

      log_prob_from_leaves[anc_index][chd_index][0] = Utils.log_sum_exp3(
                                                                        t_mat[num][0][0]+log_prob[num][0],
                                                                        t_mat[num][0][1]+log_prob[num][1],
                                                                        t_mat[num][0][2]+log_prob[num][2])

      log_prob_from_leaves[anc_index][chd_index][1] = Utils.log_sum_exp3(
                                                                        t_mat[num][1][0]+log_prob[num][0],
                                                                        t_mat[num][1][1]+log_prob[num][1],
                                                                        t_mat[num][1][2]+log_prob[num][2])

      log_prob_from_leaves[anc_index][chd_index][2] = Utils.log_sum_exp3(
                                                                        t_mat[num][2][0]+log_prob[num][0],
                                                                        t_mat[num][2][1]+log_prob[num][1],
                                                                        t_mat[num][2][2]+log_prob[num][2])
    }

    log_prob[n_index][0] = log_prob_from_leaves[n_index][0][0] + log_prob_from_leaves[n_index][1][0]
    log_prob[n_index][1] = log_prob_from_leaves[n_index][0][1] + log_prob_from_leaves[n_index][1][1]
    log_prob[n_index][2] = log_prob_from_leaves[n_index][0][2] + log_prob_from_leaves[n_index][1][2]

    p log_prob[n_index]

    log_backward[n_index][0] = -1*MAX
    log_backward[n_index][1] = 0
    log_backward[n_index][2] = -1*MAX

    # enforcely determine gg state and sg state
    if n_index == s_index
      log_backward[n_index][1] = -1*MAX
      log_backward[n_index][2] = 0
    else
      log_backward[s_index][0] = -1*MAX
      log_backward[s_index][1] = -1*MAX
      log_backward[s_index][2] = 0
    end

    last_index = @pt.sorted_nodes.size-1
    last_index.downto(0) {|num|

      if num <= @pt.n_s-1
        break
      end

      node = @pt.sorted_nodes[num]

      if gg != node && !@pt.tree.descendents(gg, root=@pt.root).include?(node)
        next
      #elsif gg == node || @pt.tree.descendents(gg, root=@pt.root).include?(node)
      #  if sg != node && !@pt.tree.descendents(sg, root=@pt.root).include?(node)
      #    next
      #  end
      end

      if node != gg
        log_backward[num][0] = Utils.log_sum_exp3(
                                                  t_mat[num][0][0]+log_backward_from_parent[num][0],
                                                  t_mat[num][1][0]+log_backward_from_parent[num][1],
                                                  t_mat[num][2][0]+log_backward_from_parent[num][2],
                                                  )
        log_backward[num][1] = Utils.log_sum_exp3(
                                                  t_mat[num][0][1]+log_backward_from_parent[num][0],
                                                  t_mat[num][1][1]+log_backward_from_parent[num][1],
                                                  t_mat[num][2][2]+log_backward_from_parent[num][2],
                                                  )
        log_backward[num][2] = Utils.log_sum_exp3(
                                                  t_mat[num][0][2]+log_backward_from_parent[num][0],
                                                  t_mat[num][1][2]+log_backward_from_parent[num][1],
                                                  t_mat[num][2][2]+log_backward_from_parent[num][2],
                                                  )

      end

      l1 = log_backward[num][0] + log_prob[num][0]
      l2 = log_backward[num][1] + log_prob[num][1]
      l3 = log_backward[num][2] + log_prob[num][2]

      state = Utils.sample_three_probs(l1, l2, l3)

      @hidden_samp[i][num] = state

      log_prob_near_leaves[num][0]     = -1*MAX
      log_prob_near_leaves[num][1]     = -1*MAX
      log_prob_near_leaves[num][2]     = -1*MAX
      log_prob_near_leaves[num][state] = 0

      @pt.tree.children(node, root=@pt.root).each do |child|
        if child == @pt.tree.root
          next
        end

        chd_index = @pt.node2num[child]

        log_backward_from_parent[chd_index][0] = log_prob_near_leaves[num][0]
        log_backward_from_parent[chd_index][1] = log_prob_near_leaves[num][1]
        log_backward_from_parent[chd_index][2] = log_prob_near_leaves[num][2]
      end

    }

    #puts "#{i} #{t_m[0][0]}"
    #p @hidden_samp
    #p log_prob

    @pt.sorted_nodes.reverse.each_with_index { |node, num|

      if num > @pt.n_s-1
        next
      end
      
      if gg != node && !@pt.tree.descendents(gg, root=@pt.root).include?(node)
        next
      elsif gg == node || @pt.tree.descendents(gg, root=@pt.root).include?(node)
        if sg != node && !@pt.tree.descendents(sg, root=@pt.root).include?(node)
          next
        end
      end

      log_backward[num][0] = Utils.log_sum_exp3(
                                                t_m[0][0]+log_backward_from_parent[num][0],
                                                t_m[1][0]+log_backward_from_parent[num][1],
                                                t_m[2][0]+log_backward_from_parent[num][2],
                                                )
      log_backward[num][1] = Utils.log_sum_exp3(
                                                t_m[0][1]+log_backward_from_parent[num][0],
                                                t_m[1][1]+log_backward_from_parent[num][1],
                                                t_m[2][2]+log_backward_from_parent[num][2],
                                                )
      log_backward[num][2] = Utils.log_sum_exp3(
                                                t_m[0][2]+log_backward_from_parent[num][0],
                                                t_m[1][2]+log_backward_from_parent[num][1],
                                                t_m[2][2]+log_backward_from_parent[num][2],
                                                )

      l1 = log_backward[num][0] + log_prob[num][0]
      l2 = log_backward[num][1] + log_prob[num][1]
      l3 = log_backward[num][2] + log_prob[num][2]

      state = Utils.sample_three_probs(l1, l2, l3)
      
      @hidden_samp[i][num] = state

      log_prob_near_leaves[num][0]     = -1*MAX
      log_prob_near_leaves[num][1]     = -1*MAX
      log_prob_near_leaves[num][2]     = -1*MAX
      log_prob_near_leaves[num][state] = 0
    }

  end

  def calc_log_Xg(g = gene, n = node)

    out_of_clade = 0

    # memo for me-> Ruby assings same reference with such declarations:
    #   Array.new(N, Array.new(K,0)). In this case, N arrays have the same reference.
    #   Confirm not to declare with above style in ruby when using multi-D array.
    log_prob_from_leaves = Array.new(@N).map{ Array.new(2).map{ Array.new(2, 0) }}
    log_prob_near_leaves = Array.new(@N).map{ Array.new(2,0) }
    log_prob             = Array.new(@N).map{ Array.new(2,0) }

    eps = 0.01
    prob_pred_error = [
                       [log(1-eps), log(eps)],
                       [log(eps), log(1-eps)]
                      ]

    # quantitate fluctuation of observed states by pre-calculated performance.
    (0..@pt.n_s-1).each{ |num|

      node  = @pt.sorted_nodes[num]
      state = @pp.profiles[@pp.symbol2index[g]][num] > 0 ? 1 : 0

      if n != node && !@pt.tree.descendents(n, root=@pt.root).include?(node)
        out_of_clade += prob_pred_error[0][state]
        next
      end

      log_prob_near_leaves[num][0] = prob_pred_error[0][state]
      log_prob_near_leaves[num][1] = prob_pred_error[1][state]
    }

    n_index = @pt.sorted_nodes.index(n)
    @pt.sorted_nodes.each_with_index { |node, num|
      if num >= n_index
        break
      end

      if n != node && !@pt.tree.descendents(n, root=@pt.root).include?(node)
        next
      end

      log_prob[num][0] = log_prob_near_leaves[num][0] + log_prob_from_leaves[num][0][0] + log_prob_from_leaves[num][1][0]
      log_prob[num][1] = log_prob_near_leaves[num][1] + log_prob_from_leaves[num][0][1] + log_prob_from_leaves[num][1][1]

      anc_index = @pt.node2num[@pt.tree.parent(node, root=@pt.root)]
      chd_index = @pt.which_child[node]

      log_prob_from_leaves[anc_index][chd_index][0] = Utils.log_sum_exp(
                                                                        @null_m[0][0]+log_prob[num][0],
                                                                        @null_m[0][1]+log_prob[num][1], false)

      log_prob_from_leaves[anc_index][chd_index][1] = Utils.log_sum_exp(
                                                                        @null_m[1][0]+log_prob[num][0],
                                                                        @null_m[1][1]+log_prob[num][1], false)


    }


    log_prob[n_index][0] = log_prob_from_leaves[n_index][0][0] + log_prob_from_leaves[n_index][1][0]
    log_prob[n_index][1] = log_prob_from_leaves[n_index][0][1] + log_prob_from_leaves[n_index][1][1]

    return log_prob[n_index][1] + out_of_clade
  end # end of method

end
