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
