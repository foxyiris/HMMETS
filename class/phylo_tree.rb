class PhyloTree

  attr_reader :sorted_nodes
  attr_reader :same_sub_tree
  attr_reader :which_child
  attr_reader :node2num
  attr_reader :tree
  attr_reader :root
  attr_reader :n_s

  def initialize(tree)

    # sorted_nodes contains leaves nodes first, then pushing ancestral nodes
    @sorted_nodes  = Array.new(0)
    @node2num      = Hash.new()
    @which_child   = Hash.new()
    @same_sub_tree = Hash.new {|h,k| h[k] = Hash.new(nil) }
    @tree          = tree
    @index         = 0
    @n_s           = 0

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
    _make_same_sub_tree(@root)

    #_add_children(root_node)

  end

  def _make_same_sub_tree(r = root)
    array = @tree.children(r, root=@root)

    if array.size != 0
      array.each do |child|
        _make_same_sub_tree(child)
        _make_same_sub_tree(child)
      end
    end

    @tree.descendents(r, root=@root).each do |desc|
      same_sub_tree[r][desc] = 1
    end

    same_sub_tree[r][r] = 1

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

=begin
# test code
if __FILE__ == $0
  require 'bundler'
  Bundler.require

  treeio = Bio::FlatFile.open(Bio::Newick, ARGV.shift)

  if newick = treeio.next_entry
    newick.options[:bootstrap_style] = :disabled
    tree = newick.tree

    pt = PhyloTree.new(tree)
  end

  naive_array = Array.new(pt.sorted_nodes.size).map{ Array.new(pt.sorted_nodes.size) }
  hash_array  = Array.new(pt.sorted_nodes.size).map{ Array.new(pt.sorted_nodes.size) }

  pt.sorted_nodes.each_with_index do |node, i|
    pt.sorted_nodes.each_with_index do |node_inside, j|
      if node == node_inside || pt.tree.descendents(node_inside, root=pt.root).include?(node)
        naive_array[i][j] = 1
      else
        naive_array[i][j] = 0
      end

      if pt.same_sub_tree[node_inside][node]
        hash_array[i][j] = 1
      else
        hash_array[i][j] = 0
      end
    end
  end

  if naive_array == hash_array
    puts "passed"
  else
    puts "not passed"
  end

end
=end
