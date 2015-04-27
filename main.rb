require 'bundler'
require_relative './class/utils.rb'
require_relative './class/phylo_prof.rb'
require_relative './class/phylo_tree.rb'
require_relative './class/infer_hist.rb'
Bundler.require

class MakePhylogenyPDF < Prawn::Document

  def initialize(pt = phylotree, gene_name, gene_branch, signal_gene_branch, hidden_states, branch_params, pp)
    super(:page_layout => :landscape)

    #gb.shift
    @pt              = pt
    @gn              = gene_name
    @gb              = gene_branch
    @sgb             = signal_gene_branch
    @hs              = hidden_states
    @bp              = branch_params
    @phylo_prof      = pp
    @root            = @pt.root
    @line_width      = 1
    @title_font_size = 20
    @font_ratio      = 1.5
    @taxon_font_size = @title_font_size/@font_ratio

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
          fill_and_stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*@font_ratio, @taxon_font_size*@font_ratio
        elsif @phylo_prof[leaf.name] == 2
          fill_color @stroke_color[3]
          stroke_color @base_stroke_color
          fill_and_stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*@font_ratio, @taxon_font_size*@font_ratio
        elsif @phylo_prof[leaf.name] == 3
          fill_color @stroke_color[2]
          stroke_color @stroke_color[2]
          fill_and_stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*@font_ratio/2, @taxon_font_size*@font_ratio

          fill_color @stroke_color[3]
          stroke_color @stroke_color[3]
          fill_and_stroke_rectangle [@leaf_offset+@prof_offset+@taxon_font_size*@font_ratio/2, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*@font_ratio/2, @taxon_font_size*@font_ratio

          stroke_color @base_stroke_color
          stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*@font_ratio, @taxon_font_size*@font_ratio
        else
          stroke_rectangle [@leaf_offset+@prof_offset, @node_pos[leaf][1]+@taxon_font_size], @taxon_font_size*@font_ratio, @taxon_font_size*@font_ratio
        end
      end
    end
  end

  def _line(r, state_parent)
    _children = @pt.tree.children(r)

    if r == @pt.root
      #@node_pos[r][0] -= 100
    end

    if _children.size != 0

      _children.each do |child|

        state = nil

        num = @pt.node2num[child]

        if !@pt.same_sub_tree[@gb[@gn]][child]
          stroke_color @stroke_color[0]
        else
          if state_parent == 0
            if child == @sgb[@gn]
              stroke_color @stroke_color[3]
              state = 2
            elsif child == @gb[@gn]
              stroke_color @stroke_color[2]
              state = 1
            else
              stroke_color @stroke_color[1]
              state = state_parent
            end
          elsif state_parent == 1

            if child == @sgb[@gn]
              stroke_color @stroke_color[3]
              state = 2
            else
              if @bp[@gn][num][0] > 0.5
                p10_hex = (255*@bp[@gn][num][0]).to_i.to_s(16)
                p10_hex = p10_hex << "0000"
                stroke_color p10_hex
                state = 0
              else
                stroke_color @stroke_color[2]
                state = state_parent
              end
            end

          elsif state_parent == 2
            if @bp[@gn][num][1] > 0.33
              p20_hex = (255*@bp[@gn][num][1]).to_i.to_s(16)
              p20_hex = p20_hex << "0000"
              stroke_color p20_hex
              state = 0
            elsif @bp[@gn][num][2] > 0.33
              p21_hex = (255*@bp[@gn][num][2]).to_i.to_s(16)
              p21_hex = "0000" << p21_hex
              stroke_color p21_hex
              state = 1
            else
              stroke_color @stroke_color[3]
              state = 2
            end
          else
            if @sgb[@gn] == child
              stroke_color @stroke_color[3]
              state = 2
            elsif @gb[@gn] == child
              stroke_color @stroke_color[2]
              state = 1
            else
              STDERR.puts "Error: unknown state is specified. #{state_parent}"
            end
          end
        end

        #if _is_signal_gained?(child)
        #  stroke_color @stroke_color[3]
        #elsif _is_gained?(child)
        #  stroke_color @stroke_color[2]
        #else
        #  stroke_color @stroke_color[0]
        #end

        stroke { line [@node_pos[r][0],@node_pos[r][1]], [@node_pos[child][0],@node_pos[child][1]] }
        _line(child, state)

        stroke_color @base_stroke_color

        draw_text "#{@pt.node2num[r]}", :at => [@node_pos[r][0],@node_pos[r][1]-@taxon_font_size/2], :size => @taxon_font_size

      end

    else
      draw_text "#{@pt.node2num[r]}", :at => [@node_pos[r][0]-@taxon_font_size*2, @node_pos[r][1]-@taxon_font_size/2], :size => @taxon_font_size
      draw_text r.name, :at => [@node_pos[r][0],@node_pos[r][1]-@taxon_font_size/2], :size => @taxon_font_size
    end

  end

  def _draw()
    max = 0
    @pt.tree.leaves(node=nil, root=@root).each do |leaf|
      d   = @pt.tree.distance(leaf, @root)
      max = d if(max < d)
    end

    stroke_color @base_stroke_color
    #stroke_axis

    text "Estimation of presequence evolution with MitoFates", :size => @title_font_size, :color=>"009900"
    draw_text "Phylogenetic profile", :at => [@leaf_offset+@prof_offset, cursor-@title_font_size]

    current_y_pos = cursor

    span = _calc_span(current_y_pos)/8

    current_y_pos -= @y_offset

    @pt.tree.leaves(node=nil, root=@root).each do |leaf|
      #draw_text "#{leaf.name}", :size => @taxon_font_size, :at => [@leaf_offset, current_y_pos]

      @node_pos[leaf] = [@leaf_offset, current_y_pos]
      current_y_pos -= span + @taxon_font_size
    end

    state = nil
    if _is_signal_gained?(@root)
      state = 2
    elsif _is_gained?(@root)
      state = 1
    else
      state = 0
    end

    _assign_node_pos(@root, max)
    _line(@root, state)
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

  # debug data
  #input_genes = ["DUMMY"]
  #input_genes = ["YDR079W"]

  # MTS data
  #input_genes = ["YNL306W"]
  #input_genes = ["YDL069C"]
  #input_genes = ["YDL069C","YPR011C"]
  input_genes = pp.symbol2index.keys
  ih = InferHist.new(pt, pp, input_genes)
  ih.pre_process()

  #ih.mcmc(30000,3000)
  ih.mcmc(30,6)

  input_genes.each do |gene|
    taxon2prof = Hash.new()
    pt.tree.leaves.each do |leaf|
      if leaf.name
        name = leaf.name.gsub(' ', '_')
        taxon2prof[leaf.name] = pp.profiles[pp.symbol2index[gene]][pp.taxon2index[name]]
      end
    end
  
    mp = MakePhylogenyPDF.new(pt, gene, ih.gain_branch, ih.signal_gain_branch, ih.hidden_states, ih.branch_params, taxon2prof)

    #mp.render_file("figs_yeast_primitive_mts/#{gene}.pdf")
    #mp.render_file("figs/#{gene}.pdf")
    #mp.render_file("figs_yeast_rbh/#{gene}.pdf")
    mp.render_file("test/#{gene}.pdf")
    #mp.render_file("#{gene}.pdf")
  end

  exit 0

end
