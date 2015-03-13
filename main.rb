require 'bundler'
Bundler.require

class Utils

  def self.log_sum_exp(a, b, f = boolean)
    return b if f
    return a + Math.log(2) if a == b

    if( a > b )
      return a + Math.log( 1 + Math.exp(b-a) )
    else
      return b + Math.log( 1 + Math.exp(a-b) )
    end

  end

  def self.log_sum_exp3(a, b, c)
    return a + Math.log(3) if a == b && b == c

    if( a > b  && a > c)
      return a + Math.log( 1 + Math.exp(b-a) + Math.exp(c-a))
    elsif ( b > a && b > c)
      return b + Math.log( 1 + Math.exp(a-b) + Math.exp(c-b))
    else
      return c + Math.log( 1 + Math.exp(a-c) + Math.exp(b-c))
    end

  end


end

class PhyloProf

  attr_reader :taxon2index
  attr_reader :profiles
  attr_reader :symbol2index

  def initialize(fp = filepath)
    @ids           = Array.new()
    @symbols       = Array.new()
    @symbol2index  = Hash.new()
    @taxons        = Array.new()
    @taxon2index   = Hash.new()
    @profiles      = Array.new{ Array.new }

    _read(fp)
  end

  def _read(fp)
    open(fp, "r"){ |file|
      index = 0
      while line = file.gets()
        line.chomp!
        if(/\AEntrez/ =~ line)
          @taxons = line.split("\t")
          @taxons.slice!(0..1)
          @taxons.each_with_index do |item, i|
            @taxon2index[item] = i.to_i
          end
        else
          @profiles[index] = line.split("\t")
          id     = @profiles[index].shift
          symbol = @profiles[index].shift
          @ids.push(id)
          @symbols.push(symbol)
          @symbol2index[symbol] = index
          index += 1
        end
      end
    }
  end

  def swap_rows(pt = phylo_tree)
    # maybe due to the reference, initialize and accessing with k and j broke array..
    # generates temp and copy its reference after each iteration is required..
    #_profiles = Array.new(@profiles.size, Array.new(@taxon2index.size, 0) )
    _profiles = Array.new()

    @profiles.each_with_index do |array, k|

      temp = Array.new()

      pt.sorted_nodes.each_with_index do |node,i|
        if i >= pt.n_s
          break
        end

        j = @taxon2index[node.name.gsub(' ', '_').to_s]
        if array[j] == "0"
          temp[i] = 0
        elsif array[j] == '1'
          temp[i] = 1
        elsif array[j] == '2'
          temp[i] = 2
        else
          STDERR.puts "Error: unexpected number for state #{array[j]}"
        end
      end
      _profiles[k] = temp
    end
    @profiles = _profiles

    pt.sorted_nodes.each_with_index do |node,i|
      @taxon2index[node.name.gsub(' ', '_').to_s] = i
    end

  end
  
end


class PhyloTree

  attr_reader :sorted_nodes
  attr_reader :which_child
  attr_reader :node2num
  attr_reader :tree
  attr_reader :root
  attr_reader :n_s

  def initialize(tree)

    # sorted_nodes contains leaves nodes first, then pushing ancestral nodes
    @sorted_nodes = Array.new(0)
    @node2num     = Hash.new()
    @which_child  = Hash.new()
    @tree         = tree
    @index        = 0
    @n_s          = 0

    # leaves returns root of root due to the link number, so omit it.
    # in other words, #leaf != tree.leaves.size
    @tree.leaves.each do |leaf|
      if leaf == @tree.root
        next
      else
        @n_s += 1
      end
    end
    
    @leaf_index = 0
    @int_index  = @n_s

    # Caution: tree.root specifies ancestor of root.
    # so children(root) returns true root.
    # In fact, @tree.descendents(@tree.root).size == 2S-1
    @root = @tree.children(@tree.root)[0]

    _make_sorted_array(@root)
    _make_hash(@sorted_nodes)

    #_add_children(root_node)

  end

  # make an array which contains leaves in first @n_s pos,
  # and internal nodes in the left @n_s-1 positions
  # plus, subnode contains smaller number so parsing is possible from bottom.
  def _make_sorted_array(node)
    flag = Hash.new()

    array = @tree.children(node)
    
    if array.size != 0
      i = 0
      array.each do |child|
        _make_sorted_array(child)
        @which_child[child] = i
        i += 1
      end
    end

    if array.size == 0
      #STDERR.puts "leaf: #{@leaf_index} #{node.name}"
      @sorted_nodes[@leaf_index] = node
      @leaf_index += 1
    else
      #STDERR.puts "inte: #{@int_index} #{node.name}"
      @sorted_nodes[@int_index] = node
      @int_index += 1
    end

  end

  def _make_hash(array)
    if array.empty?
      exit 1
    end

    array.each_with_index do |item, i|
      @node2num[item] = i
    end

  end

  # deprecated
  def _add_children(node)
    begin
      array = @tree.children(node)
    rescue Exception => e
      p e.message
    end

    if array.size != 0
      # internal node
      _add_children(array[0])
      _add_children(array[1])
      @heap[@index+=1] = node
    else
      # leave
      @heap[@index+=1] = node
    end
  end

end


#FP_rate = 0.005
#FN_rate = 0.005

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
    @is_clade    = Array.new(@N).map{ Array.new(@N) }
    @null_m      = [
                    [log(1), log(0)],
                    [log(0.03), log(0.97)]
                   ]


    # in the future this has to be k-D array? to store MULTI gainings.
    @gain_branch        = Hash.new()
    @signal_gain_branch = Hash.new()
    # lost branch can be many, so they are contained as an array
    @gene_lost_branch   = Hash.new {|hash,key| hash[key] = Array.new()}
    @signal_lost_branch = Hash.new {|hash,key| hash[key] = Array.new()}

    # this is an array which stores the names of the target genes
    @input_genes  = input_genes
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

    p most_likel_state_vec

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

class MakePhylogenyPDF < Prawn::Document

  def initialize(pt = phylotree, gene_name, gene_branch, signal_gene_branch, pp)
    super(:page_layout => :landscape)

    #gb.shift
    @pt              = pt
    @gn              = gene_name
    @gb              = gene_branch
    @sgb             = signal_gene_branch
    @phylo_prof      = pp
    @root            = @pt.root
    @line_width      = 1
    @title_font_size = 20
    @taxon_font_size = @title_font_size/2

    @leaf_offset = 350
    @y_offset    = 50
    @x_offset    = 50
    @prof_offset = 150

    @node_pos    = Hash.new{|hash, key| hash[key] = Array.new(2)}

    @base_stroke_color = "333333"
    @stroke_color = ["CCCCCC", "ff0000", "0000ff", "00ff00"] # out of gain, loss, gene gain and signal gain

    line_width @line_width
    #font_size @font_size
    _draw()
  end

  def _calc_span(p = current_y_pos)
    span = 10 # dummy number to avoid crush

    weight=1
    while(span < 11)
      @taxon_font_size *= weight
      range = p - 2*@y_offset
      range -= @pt.n_s * (@taxon_font_size * weight)     
      span = range / (@pt.n_s-1)
      weight -= 0.1
    end

    return span

  end

  def _is_gained?(node)
    return @pt.tree.descendents(@gb[@gn], root=@root).include?(node) || @gb[@gn] == node
  end

  def _is_signal_gained?(node)
    return @pt.tree.descendents(@sgb[@gn], root=@root).include?(node) || @sgb[@gn] == node
  end

  # recursion
  # x position is determined by max and defined range, say 350
  # y position is determined by leaf position.
  def _assign_node_pos(r, max)
    _children = @pt.tree.children(r, root=@root)
    unit = (@leaf_offset-@x_offset)/max

    if _children.size != 0
      #internal node
      x_pos  = 0

      dist = @pt.tree.distance(r, @root)
      
      x_pos  = @x_offset + dist*unit      

      if @node_pos.key?(_children[0]) && !@node_pos.key?(_children[1])
        _assign_node_pos(_children[1], max)
      elsif !@node_pos.key?(_children[0]) && @node_pos.key?(_children[1])
        _assign_node_pos(_children[0], max)
      else
        _assign_node_pos(_children[0], max)
        _assign_node_pos(_children[1], max)
      end

      if @node_pos.key?(_children[0]) && @node_pos.key?(_children[1])
        y_pos  = 0

        _children.each do |child|
          y_pos += @node_pos[child][1]
          #x_dist.push(@pt.tree.distance(child, @pt.tree.parent(child)))
        end

        @node_pos[r] = [x_pos, y_pos/2]

      end

    else
      # leaf
      if !@node_pos.key?(r)
        STDERR.puts("Error: #{r.name} has no position.")
      end
    end

  end

  def _add_root_annotation
    if _is_signal_gained?(@root)
      stroke_color @stroke_color[3]
    elsif _is_gained?(@root)
      stroke_color @stroke_color[2]
    else
      stroke_color @stroke_color[0]
    end

    stroke{ line [@node_pos[@root][0]-15, @node_pos[@root][1]],@node_pos[@root]}

    _width = width_of("Last Common Ancestor", :size => @taxon_size_font)
    text_box("Last Common Ancestor",
             :at => [@node_pos[@root][0]-30, @node_pos[@root][1] - _width/2],
             :rotate => 90, :rotate_around => [@node_pos[@root][0]-30, @node_pos[@root][1] - _width/2])

    stroke_color @base_stroke_color
  end

  def _draw_prof
    @pt.tree.leaves(node=nil, root=@root).each do |leaf|
      if leaf.name
        if @phylo_prof[leaf.name] == 1
          fill_color @stroke_color[2]
          stroke_color @base_stroke_color
          fill_and_stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*2, @taxon_font_size*2
        elsif @phylo_prof[leaf.name] == 2
          fill_color @stroke_color[3]
          stroke_color @base_stroke_color
          fill_and_stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*2, @taxon_font_size*2
        else
          stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*2, @taxon_font_size*2
        end
      end
    end
  end

  def _line(r)
    _children = @pt.tree.children(r)

    if r == @pt.root
      #@node_pos[r][0] -= 100
    end

    if _children.size != 0

      _children.each do |child|
        if _is_signal_gained?(child)
          stroke_color @stroke_color[3]
        elsif _is_gained?(child)
          stroke_color @stroke_color[2]
        else
          stroke_color @stroke_color[0]
        end

        stroke { line [@node_pos[r][0],@node_pos[r][1]], [@node_pos[child][0],@node_pos[child][1]] }
        _line(child)

        stroke_color @base_stroke_color

      end

    else
      draw_text r.name, :at => [@node_pos[r][0],@node_pos[r][1]-@taxon_font_size/2]
    end

  end

  def _draw()
    max = 0
    @pt.tree.leaves(node=nil, root=@root).each do |leaf|
      d   = @pt.tree.distance(leaf, @root)
      max = d if(max < d)
    end

    stroke_color @base_stroke_color
    stroke_axis

    text "Estimation of presequence evolution with MitoFates", :size => @title_font_size, :color=>"009900"
    draw_text "Phylogenetic profile", :at => [@leaf_offset+@prof_offset, cursor-@title_font_size]

    current_y_pos = cursor

    span = _calc_span(current_y_pos)

    current_y_pos -= @y_offset

    @pt.tree.leaves(node=nil, root=@root).each do |leaf|
      #draw_text "#{leaf.name}", :size => @taxon_font_size, :at => [@leaf_offset, current_y_pos]

      @node_pos[leaf] = [@leaf_offset, current_y_pos]
      current_y_pos -= span + @taxon_font_size
    end

    _assign_node_pos(@root, max)
    _line(@root)
    _add_root_annotation
    _draw_prof
    

  end

end

# main 
if __FILE__ == $0
=begin
  require 'optparse'
  Version = "0.0.2"
  args = {}
  OptionParser.new do |parser|
    parser.on('-u', '--unikont', 'Flag whether amoeba is treated as unikont or independent category amoeba.'){|v| args[:unikont] = v}
    parser.permute!(ARGV)
  end
=end


  treeio = Bio::FlatFile.open(Bio::Newick, ARGV.shift)
  pp     = PhyloProf.new(ARGV.shift)
  pt     = nil

  if newick = treeio.next_entry
    newick.options[:bootstrap_style] = :disabled
    tree = newick.tree

    pt = PhyloTree.new(tree)
  end

  pp.swap_rows(pt)

  input_genes = ["YDL069C"]
  #input_genes = pp.symbol2index.keys
  ih = InferHist.new(pt, pp, input_genes)
  ih.pre_process()

  #if ih.gain_branch[input_genes[0]].name != ""
  #  p ih.gain_branch[input_genes[0]].name
  #else
  #  p pt.tree.descendents(ih.gain_branch[input_genes[0]])
  #end
  #p pp.profiles[pp.symbol2index[input_genes[0]]]

  taxon2prof = Hash.new()
  pt.tree.leaves.each do |leaf|
    if leaf.name
      name = leaf.name.gsub(' ', '_')
      taxon2prof[leaf.name] = pp.profiles[pp.symbol2index[input_genes[0]]][pp.taxon2index[name]]
    end
  end
  
  mp = MakePhylogenyPDF.new( pt, input_genes[0], ih.gain_branch, ih.signal_gain_branch, taxon2prof)
  mp.render_file("prawn.pdf")

  exit 0
end
