require 'fileutils'
require 'bundler'
require_relative './class/utils.rb'
require_relative './class/phylo_prof.rb'
require_relative './class/phylo_tree.rb'
require_relative './class/infer_hist.rb'
Bundler.require

# main 
if __FILE__ == $0

  use_multicore = true
  if use_multicore
    processor_count = Parallel.processor_count
  else
    processor_count = 1
  end

  model_file_path = "/home/yoshinori/PhyloFates/models/"

  if !Dir.exists?(model_file_path)
    STDERR.puts "Error: #{model_file_path} doesn't exits. Make it."
    FileUtils.mkdir_p(model_file_path)
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

  input_genes = pp.symbol2index.keys

  ih = InferHist.new(pt, pp, input_genes, processor_count)
  ih.pre_process_parallel()

  total_step = 50000
  burn_in    = 3000

  ih.mcmc(total_step, burn_in, model_file_path)

  exit 0

end
