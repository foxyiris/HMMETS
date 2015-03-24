class Utils

  def self.sample_three_probs(log1, log2, log3)
    max = [log1,log2,log3].max

    p1 = Math.exp(log1-max)
    p2 = Math.exp(log2-max)
    p3 = Math.exp(log3-max)

    sum = p1+p2+p3
    p1 /= sum
    p2 /= sum
    p3 /= sum

    prng = Random.new()
    random = prng.rand(1.0)

    if random < p1
      return 0
    elsif p1 <= random && random < p2
      return 1
    elsif p2 <= random
      return 2
    else
      return -1 #error check
    end
  end

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
    return self.log_sum_exp(self.log_sum_exp(a, b, false), c, false)
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
