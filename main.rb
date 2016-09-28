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
  #input_genes = ["Matryoshka_3692"]

  #input_genes = ["Matryoshka_1156", "Matryoshka_1578", "Matryoshka_3348", "Matryoshka_3527", "Matryoshka_1472", "Matryoshka_2999", "Matryoshka_2906", "Matryoshka_3722", "Matryoshka_1518", "Matryoshka_3152", "Matryoshka_2858", "Matryoshka_2373", "Matryoshka_2443", "Matryoshka_3432", "Matryoshka_3149", "Matryoshka_3089", "Matryoshka_1844", "Matryoshka_2094", "Matryoshka_2522", "Matryoshka_3692", "Matryoshka_3555", "Matryoshka_1533", "Matryoshka_2856", "Matryoshka_2949", "Matryoshka_3740", "Matryoshka_1262", "Matryoshka_2953", "Matryoshka_2357", "Matryoshka_1272", "Matryoshka_2793", "Matryoshka_2150", "Matryoshka_2743", "Matryoshka_2781", "Matryoshka_3662", "Matryoshka_2850", "Matryoshka_3399", "Matryoshka_1515", "Matryoshka_2449", "Matryoshka_3264", "Matryoshka_1573", "Matryoshka_3002", "Matryoshka_3530", "Matryoshka_3338", "Matryoshka_1475", "Matryoshka_1535", "Matryoshka_2225", "Matryoshka_2340", "Matryoshka_1903", "Matryoshka_3410", "Matryoshka_2476", "Matryoshka_1430", "Matryoshka_3336", "Matryoshka_2289", "Matryoshka_3628", "Matryoshka_3180", "Matryoshka_1754", "Matryoshka_2608", "Matryoshka_3562", "Matryoshka_2914", "Matryoshka_1839", "Matryoshka_3155", "Matryoshka_3637", "Matryoshka_10439", "Matryoshka_1770", "Matryoshka_2706", "Matryoshka_3276", "Matryoshka_16127", "Matryoshka_2228", "Matryoshka_1825", "Matryoshka_1848", "Matryoshka_10743", "Matryoshka_2090", "Matryoshka_2526", "Matryoshka_3455", "Matryoshka_1961", "Matryoshka_1958", "Matryoshka_1206", "Matryoshka_2602", "Matryoshka_3129", "Matryoshka_3165", "Matryoshka_15202", "Matryoshka_3658", "Matryoshka_3177", "Matryoshka_2885", "Matryoshka_1314", "Matryoshka_1809", "Matryoshka_2063", "Matryoshka_2661", "Matryoshka_2515", "Matryoshka_2853", "Matryoshka_2505", "Matryoshka_1364", "Matryoshka_2091", "Matryoshka_3206", "Matryoshka_1758", "Matryoshka_1051", "Matryoshka_2741", "Matryoshka_1715", "Matryoshka_11178", "Matryoshka_1384", "Matryoshka_2202", "Matryoshka_1099", "Matryoshka_1244", "Matryoshka_1137", "Matryoshka_3374"]
  #input_genes = ["Matryoshka_3581", "Matryoshka_3586", "Matryoshka_3591", "Matryoshka_3655", "Matryoshka_3668", "Matryoshka_3675", "Matryoshka_3742", "Matryoshka_3800", "Matryoshka_3880", "Matryoshka_3882", "Matryoshka_3901", "Matryoshka_3909", "Matryoshka_3922", "Matryoshka_3969", "Matryoshka_3972", "Matryoshka_3989", "Matryoshka_4005", "Matryoshka_4017", "Matryoshka_4050", "Matryoshka_4075", "Matryoshka_4095", "Matryoshka_4109", "Matryoshka_4123", "Matryoshka_4130", "Matryoshka_4171", "Matryoshka_4180", "Matryoshka_4181", "Matryoshka_4189", "Matryoshka_4206", "Matryoshka_4246", "Matryoshka_4251", "Matryoshka_4296", "Matryoshka_4319", "Matryoshka_4350", "Matryoshka_4380", "Matryoshka_4446", "Matryoshka_4504", "Matryoshka_4559", "Matryoshka_4578", "Matryoshka_4603", "Matryoshka_4606", "Matryoshka_4607", "Matryoshka_4682", "Matryoshka_4697", "Matryoshka_4742", "Matryoshka_4806", "Matryoshka_4842", "Matryoshka_4845", "Matryoshka_4850", "Matryoshka_4928", "Matryoshka_4974", "Matryoshka_5004", "Matryoshka_5058", "Matryoshka_5072", "Matryoshka_5179", "Matryoshka_5202", "Matryoshka_5236", "Matryoshka_5346", "Matryoshka_5418", "Matryoshka_5432", "Matryoshka_5469", "Matryoshka_5473", "Matryoshka_5490", "Matryoshka_5590", "Matryoshka_5698", "Matryoshka_5787", "Matryoshka_5983", "Matryoshka_5985", "Matryoshka_6011", "Matryoshka_6018", "Matryoshka_6133", "Matryoshka_6167", "Matryoshka_6181", "Matryoshka_6371", "Matryoshka_6451", "Matryoshka_6664", "Matryoshka_6694", "Matryoshka_6753", "Matryoshka_6754", "Matryoshka_6889", "Matryoshka_7010", "Matryoshka_7042", "Matryoshka_7270", "Matryoshka_7554", "Matryoshka_7653", "Matryoshka_7748", "Matryoshka_7757", "Matryoshka_7783", "Matryoshka_7992", "Matryoshka_8333", "Matryoshka_8472", "Matryoshka_8867", "Matryoshka_8874", "Matryoshka_8928", "Matryoshka_9835"]

  input_genes = pp.symbol2index.keys

  # insufficient in 273 Leca set
  #input_genes = ["Matryoshka_2090", "Matryoshka_2202", "Matryoshka_2225", "Matryoshka_2228", "Matryoshka_2289", "Matryoshka_2515", "Matryoshka_2602", "Matryoshka_2643", "Matryoshka_2706", "Matryoshka_2850", "Matryoshka_2906", "Matryoshka_2961", "Matryoshka_3097", "Matryoshka_3165", "Matryoshka_3177", "Matryoshka_3180", "Matryoshka_3291", "Matryoshka_3348", "Matryoshka_3527", "Matryoshka_3586", "Matryoshka_3637", "Matryoshka_3751", "Matryoshka_3791", "Matryoshka_3909", "Matryoshka_3922", "Matryoshka_3969", "Matryoshka_4043", "Matryoshka_4181", "Matryoshka_4246", "Matryoshka_4250", "Matryoshka_4335", "Matryoshka_4443", "Matryoshka_4444", "Matryoshka_4446", "Matryoshka_4504", "Matryoshka_4697", "Matryoshka_4729", "Matryoshka_4789", "Matryoshka_4799", "Matryoshka_4806", "Matryoshka_4976", "Matryoshka_5072", "Matryoshka_5179", "Matryoshka_5185", "Matryoshka_5236", "Matryoshka_5432", "Matryoshka_5518", "Matryoshka_5575", "Matryoshka_5983", "Matryoshka_6018", "Matryoshka_6074", "Matryoshka_6181", "Matryoshka_6753", "Matryoshka_6754", "Matryoshka_7010", "Matryoshka_7042", "Matryoshka_7266", "Matryoshka_7270", "Matryoshka_7992", "Matryoshka_8867", "Matryoshka_8874"]

  ih = InferHist.new(pt, pp, input_genes, processor_count)
  ih.pre_process_parallel()

  #honban
=begin
  total_step = 50000
  burn_in    = 3000
=end

  # debug
  total_step = 25000
  burn_in    = 1500

  #ih.mcmc(total_step, burn_in, "/home/yoshinori/MitoFates_MCMC/MitoFatesProb/parallel_param_files/Recalc_post_sg/")
  #ih.mcmc(100000, 3000)
  #ih.mcmc(10, 4)

  exit 0

end
