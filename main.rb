require 'bundler'
require_relative './class/utils.rb'
require_relative './class/phylo_prof.rb'
require_relative './class/phylo_tree.rb'
require_relative './class/infer_hist.rb'
Bundler.require

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

  processor_count = Parallel.processor_count
  use_multicore = true
  if use_multicore
    processor_count -= 10
  end

  treeio = Bio::FlatFile.open(Bio::Newick, ARGV.shift)
  pp     = PhyloProf.new(ARGV.shift)
  pt     = nil

  if newick = treeio.next_entry
    newick.options[:bootstrap_style] = :disabled
    tree = newick.tree

    pt = PhyloTree.new(tree)
  end

  pp.swap_rows(pt)

=begin
  if true
    pt.sorted_nodes.each_with_index do |node, i|
      anc_index = pt.node2num[pt.tree.parent(node, root=pt.root)]
      puts "#{i} #{anc_index}"
    end
    exit 0
  end
=end

  # debug data
  #input_genes = ["DUMMY"]
  #input_genes = ["YDR079W"]
  #input_genes = ["2A5D_YEAST"]

  # MTS data
  #input_genes = ["YNL306W"]
  #input_genes = ["YDL069C"]
  #input_genes = ["YDL069C","YPR011C"]
  input_genes = pp.symbol2index.keys
  ih = InferHist.new(pt, pp, input_genes, processor_count)
  ih.pre_process()

  total_step = 100000
  burn_in    = 3000

  ih.mcmc(total_step, burn_in, "/home/yoshinori/MitoFates_MCMC/MitoFatesProb/parallel_param_files/TIM50/")
  #ih.mcmc(100000, 3000)
  #ih.mcmc(10, 4)

  exit 0

end
