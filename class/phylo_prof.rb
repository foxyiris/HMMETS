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
