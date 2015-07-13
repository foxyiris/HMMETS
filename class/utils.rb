class Utils
  require 'bundler/setup'
  require 'inline'

  inline do |builder|
    builder.include('<math.h>')
    builder.c_singleton'
      double
      log_c(double i){
        return log(i);
      }'
  end

  inline do |builder|
    builder.include('<math.h>')
    builder.c_singleton'
      double
      exp_c(double i){
        return exp(i);
      }'
  end

  inline do |builder|
    builder.include('<math.h>')
    builder.c_singleton'
      double
      log_sum_exp_c(double a, double b){
        if(a == b)
          return a + log(2);

        if( a > b )
          return a + log( 1 + exp(b-a) );
        else
          return b + log( 1 + exp(a-b) );
      }'
  end

  def self.calc_mean(array)
    sum = 0.0
    array.each {|v|
      sum += v.to_f
    }
    return sum/array.size.to_f
  end

  def self.calc_var(array)
    mu = self.calc_mean(array)
    sum = 0.0
    array.each {|v|
      sum += (mu - v.to_f)**2.0
    }
    return sum/(array.size - 1).to_f
  end

  def self.calc_auto_cor_time(vals)
    # I use batch means to estimate autocorrelation time tau.
    m      = vals.length**(2/3.to_f) # size of batch
    nbatch = vals.length**(1/3.to_f) # the number of batch

    batch_means = Array.new()
    vals.each_slice(m.to_i).to_a.each do |batch|
      next if batch.size < m.to_i
      batch_means.push(self.calc_mean(batch))
    end

    sigma_m = self.calc_var(batch_means)

    if self.calc_var(vals) == 0
      return (sigma_m*m.to_i)/1e-10
    else
      return (sigma_m*m.to_i)/(self.calc_var(vals))
    end
  end

  def self.calc_ess(vals)
    tau = self.calc_auto_cor_time(vals)
    return vals.size/tau.to_f
  end

  def self.sample_three_probs(log1, log2, log3)
    max = [log1,log2,log3].max

    p1 = exp_c(log1-max)
    p2 = exp_c(log2-max)
    p3 = exp_c(log3-max)

    sum = p1+p2+p3
    p1 /= sum
    p2 /= sum
    p3 /= sum

    prng = Random.new()
    random = prng.rand(1.0)
    if random < p1
      return 0
    elsif p1 <= random && random < p1+p2
      return 1
    elsif p1+p2 <= random && random < p1+p2+p3
      return 2
    else
      return -1 #error check
    end
  end

  def self.log_sum_exp(a, b, f = boolean)
    return b if f
    return a + log_c(2) if a == b

    if( a > b )
      return a + log_c( 1 + exp_c(b-a) )
    else
      return b + log_c( 1 + exp_c(a-b) )
    end

  end

  def self.log_sum_exp3(a, b, c)
    return self.log_sum_exp_c(self.log_sum_exp_c(a, b), c)
  end

=begin
  # deprecated. makes NaN error sometimes
  def self.log_sum_exp3(a, b, c)
    return a + Math.log(3) if a == b && b == c

    if( a > b  && a > c)
      v = a + Math.log( 1 + Math.exp(b-a) + Math.exp(c-a))
      return v
    elsif ( b > a && b > c)
      v = b + Math.log( 1 + Math.exp(a-b) + Math.exp(c-b))
      return v
    else
      v = c + Math.log( 1 + Math.exp(a-c) + Math.exp(b-c))
      return v
    end

  end
=end

end

=begin
# test code 
if __FILE__ == $0
  puts Utils.log_c()
  puts Math.log(0.2)
end

# test for sampling
if __FILE__ == $0
  100.times do
    puts Utils.sample_three_probs(-5.415927617733775, -4.724336573224139, Utils.log_c(0))
  end
end
=end

# test ess
=begin
if __FILE__ == $0
  # 29.7017
  vars =   "-0.4086230949 -1.3397993775  0.0335816352 -0.0578464063  1.0518453543 -0.3547756144 -2.0671772765  1.0159316828  1.1376771106  0.8089929164 -0.4693850321 0.6488526393 -0.0046054833  0.4671574505  0.3315279091  0.2393821731  0.3210750351 -1.5675064806 -1.1262958177  0.1031456429 -0.2096699395 -0.4903965588 0.6788403590  0.8200379558  0.2195079622  1.1838240212 -0.3259121181  0.2763020061  0.4492294543 -0.7190543794 -0.7839148379  1.7118517262  0.9106177494 -1.1267160890  0.0226193251  0.8636777040  0.5385234803  1.3219670511  1.0334817805  1.5006690972 -0.6783059738 -0.3742272096  1.3127189882  0.3373372029 2.4428805754  0.0008477733  0.0246747578  0.2209306756  0.9277670570 -2.0448757971  0.8988645250  1.6996891669  1.0939219288 -0.1723412777 -0.3386250686 -0.1983003995  1.0105133580  0.3161174880 -0.7070337209  0.2248077001  0.8703202908  1.5961459270  1.8763018723 -1.1936848829 -0.0976222139  0.5679758301 -0.2849732565 -1.4121658578  1.5012210202  0.3459440341 -0.2227936064 -1.2157065874 -1.3268707195  0.0430743464 -1.1666242576 -0.8807011724 -0.1954612072 0.1718725545  0.4007181720 -1.1634983165  1.1073717318  0.1863529393 -0.2327152917 -1.0868135077  0.5931523549  0.7320514190  0.1931042244  0.0915942367 -0.6904350982 -1.6261701467  0.3671645283 -0.3337950454 -1.0277899944 -0.5187557401  2.2118148559  2.1867467842 -0.1527743645  0.3337521042  1.2223470856 0.5978201594".split(" ")
  p Utils.calc_ess(vars)
end
=end
