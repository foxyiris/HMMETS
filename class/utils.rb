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
