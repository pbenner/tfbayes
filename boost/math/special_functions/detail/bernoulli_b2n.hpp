
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2013 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_BERNOULLI_B2N_2013_05_30_HPP_
 #define _BOOST_BERNOULLI_B2N_2013_05_30_HPP_

  #include <limits>
  #include <vector>
  #include <cmath>
  #include <boost/cstdint.hpp>
  #include <boost/math/special_functions/pow.hpp>
  #include <boost/math/policies/policy.hpp>
  #include <boost/math/policies/error_handling.hpp>
  #include <boost/static_assert.hpp>
  #include <boost/mpl/if.hpp>
  #include <boost/mpl/int.hpp>
  #include <boost/type_traits/is_convertible.hpp>
  #include <boost/multiprecision/cpp_dec_float.hpp>
  #include <boost/math/constants/constants.hpp>

  using namespace boost::multiprecision;
  using std::size_t;

  struct pol{};

  namespace boost { namespace math { namespace detail {

  template<class T>
  bool bernouli_impl_index_does_overflow(const int n)
  {
    // There are certain cases for which the index n, when trying
    // to compute Bn, is known to overflow. In particular, when
    // using 32-bit float and 64-bit double (IEEE 754 conformant),
    // overflow will occur if the index exceeds the amount allowed
    // in the tables of Bn.

    bool the_index_does_overflow = false;

    if(std::numeric_limits<T>::is_specialized)
    {
      if(std::numeric_limits<T>::max_exponent == 128)
      {
        // This corresponds to 4-byte float, IEEE 745 conformant.
        the_index_does_overflow = (n >= 64);
      }

      if(std::numeric_limits<T>::max_exponent == 1024)
      {
        // This corresponds to 8-byte float, IEEE 745 conformant.
        the_index_does_overflow = (n >= 260);
      }
    }

    return the_index_does_overflow;
  }

  template<class T>
  bool bernouli_impl_index_might_overflow(const int n)
  {
    if(std::numeric_limits<T>::is_specialized == false)
    {
      // If numeric limits is not specialized, then not much is known
      // about overflow of Bn. So we examine the small-index behavior
      // of Bn and notice that |Bn| > 1 for n >= 14.
      return (n >= 14);
    }

    if(bernouli_impl_index_does_overflow<T>(n))
    {
      // If the index *does* overflow, then it also *might* overflow.
      return true;
    }

    // Here, we use an asymptotic expansion of |Bn| from Luschny
    // to estimate if a given index n for Bn *might* overflow.
    const T log_of_four_pi = T(log(boost::math::constants::two_pi<T>() * 2));
    const T two_pi_e       = boost::math::constants::two_pi<T>() * boost::math::constants::e<T>();
    const float nf = static_cast<float>(n);
    const T nx (nf);
    const T nx2(nx * nx);
    const T n_log_term     = (nx + boost::math::constants::half<T>()) * T(log(nx / two_pi_e));

    const T approximate_log_of_bn =   log_of_four_pi
                                    + n_log_term
                                    + boost::math::constants::half<T>()
                                    + (T(1) / (nx * 12))
                                    - (T(1) / (nx2 * nx * 360))
                                    + (T(1) / (nx2 * nx2 * nx * 1260))
                                    + n * boost::math::constants::ln_two<T>()
                                    - log(nx)
                                    + T(log(ldexp(T(1), n) - 1));

    const T approximate_exponent2_of_bn =   approximate_log_of_bn / boost::math::constants::ln_two<T>();

    static const T max_exponent2(std::numeric_limits<T>::max_exponent);

    // Check if the index might overflow. :Here we multiply the estimated
    // value once again by 1.1 in order to remain conservative with
    // the prediction of potential overflow.
    const bool the_index_might_overflow = (T(approximate_exponent2_of_bn * 1.1F) > max_exponent2);

    return the_index_might_overflow;
  }

  template<class T>
  int possible_overflow_index()
  {
    // we use binary search to determine a good approximation for an indec that might overflow

    int upper_limit = 10000;
    int lower_limit = 1;

    if(bernouli_impl_index_might_overflow<T>(upper_limit * 2) == 0)
    {
      return upper_limit;
    }

    while(upper_limit > (lower_limit + 4))
    {
      const int mid = (upper_limit + lower_limit) / 2;

      if(bernouli_impl_index_might_overflow<T>(mid * 2) == 0)
      {
        lower_limit = mid;
      }
      else
      {
        upper_limit = mid;
      }
    }

    return lower_limit;
  }

  template <class T>
  struct max_bernoulli_index;

  template<class T>
  inline T unchecked_bernoulli_b2n(size_t n);

  template<class T,class TypeIterator,class Policy>
  inline void tangent(TypeIterator tangent_numbers,const int &m, T , Policy&)
  {
    // TBD: Initialize this before main().
    static size_t min_overflow_index = static_cast<size_t>(possible_overflow_index<T>());

    tangent_numbers[0U] = T(0U);
    tangent_numbers[1U] = T(1U);

    for(size_t k = 2U; k <= size_t(m); k++)
    {
      if(   (k >= min_overflow_index)
         && (boost::math::tools::max_value<T>()/(k - 1) < tangent_numbers[k - 1])
        )
      {
        tangent_numbers[k]=policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>", "Overflow error while calculating tangent number %1%",Policy());
      }
      else
      {
        tangent_numbers[k] = (T(k - 1)) * tangent_numbers[k - 1];
      }
    }

    for(size_t k = 2; k <= size_t(m); k++)
    {
      for(size_t j = k; j <= size_t(m); j++)
      {
        if(   (j >= min_overflow_index)
           && (   (boost::math::tools::max_value<T>() / (j - k) < tangent_numbers[j - 1])
               || (boost::math::tools::max_value<T>() / (j - k + 2) < tangent_numbers[j])
               || (boost::math::tools::max_value<T>() - tangent_numbers[j - 1] * (j - k) < tangent_numbers[j] * (j - k + 2))
               || ((boost::math::isinf)(tangent_numbers[j])))
          )
        {
          tangent_numbers[j] = policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>", "Overflow error while calculating tangent number %1%",Policy());
        }
        else
        {
          tangent_numbers[j] = (tangent_numbers[j - 1] * (j - k)) + (tangent_numbers[j] * (j - k + 2));
        }
      }
    }
  }

  template <class T, class Policy>
  T tangent_numbers(const int nn, Policy &pol)
  {
    if((nn % 2) != 0)
    {
      T x(0U);
      return x;
    }

    const boost::int32_t m = nn / 2;

    std::vector<T> tangent_numbers(m + 1);

    boost::math::detail::tangent(tangent_numbers.begin(),m,T(1),pol);

    T power_two(1);
    power_two = ldexp(T(1), 2 * m);

    T x(nn);
    x  = x / (power_two*(power_two-1));
    x *= tangent_numbers[nn / 2];

    return (((nn % 4) == 0) ? -x : x);
  }

  template <class T, class Policy>
  void tangent_numbers_series(std::vector<T>& bn, const size_t start_index, const size_t number_of_bernoullis_bn, Policy &pol)
  {
    const size_t m = start_index + number_of_bernoullis_bn;

    std::vector<T> tangent_numbers(m + 1U);

    const int mn = static_cast<int>(m);
    boost::math::detail::tangent(tangent_numbers.begin(), mn, T(1), pol);

    T power_two(1);

    power_two = ldexp(T(1), 2 * static_cast<boost::int32_t>(start_index));

    bn.clear();
    bn.resize(number_of_bernoullis_bn);

    for(size_t i = 0; i < number_of_bernoullis_bn; i++)
    {
      T b((T(i + start_index) * 2));

      b  = b / (power_two * (power_two - 1));
      b *= tangent_numbers[i+start_index];

      power_two *= 4;

      const bool b_neg = (static_cast<boost::int32_t>((i + start_index) % static_cast<boost::int32_t>(2)) == static_cast<boost::int32_t>(0));

      bn[i] = ((!b_neg) ? b : -b);
    }
  }

  template <class T, class OutputIterator, class Policy>
  inline OutputIterator cache_imp(size_t start_index,
                                      size_t number_of_bernoullis_bn,
                                      OutputIterator out_it,
                                      const Policy& pol)
  {

    if((start_index + number_of_bernoullis_bn - 1) <= max_bernoulli_index<T>::value)
    {
      OutputIterator last= out_it + number_of_bernoullis_bn;

      while(out_it!=last)
      {
        *out_it = unchecked_bernoulli_b2n<T>(start_index);
        ++out_it;
        ++start_index;
      }

      //return one past the last element
      return out_it;
    }
    else if(   ( start_index <= static_cast<boost::int32_t>(max_bernoulli_index<T>::value))
            && ((start_index + number_of_bernoullis_bn) > static_cast<boost::int32_t>(max_bernoulli_index<T>::value))
           )
    {
      out_it=cache_imp<T,OutputIterator,Policy>(start_index,max_bernoulli_index<T>::value - start_index +1, out_it,pol);

      out_it=cache_imp<T,OutputIterator,Policy>(max_bernoulli_index<T>::value+1,
                                                number_of_bernoullis_bn - max_bernoulli_index<T>::value + start_index -1,
                                                out_it,
                                                pol);
      return out_it;
    }

    std::vector<T> bn;

    tangent_numbers_series(bn, start_index , number_of_bernoullis_bn, pol);

    OutputIterator last = out_it + number_of_bernoullis_bn;

    size_t i = 0U;

    while(out_it != last)
    {
      *out_it = bn[i];
      ++out_it;
      ++i;
    }

    return out_it;
  }

  template <class T, class Policy>
  T bernoulli_number_imp(const int n, const Policy &pol)
  {
    if(n < 0)
    {
      policies::raise_domain_error<T>("boost::math::bernoulli<%1%>", "Index should be >= 0 but got %1%", T(n/2), Policy());
    }

    static std::vector<T> cache_table;
    std::size_t previous_cache_size = cache_table.size();
    std::size_t index = n / 2;

    if( previous_cache_size <= index)
    {
      cache_table.resize(index + 1);

      cache_imp<T>(previous_cache_size,
                   index - previous_cache_size + 1,
                   cache_table.begin() + previous_cache_size,
                   pol);
    }
    return cache_table[index];
  }

  template <class T, class OutputIterator, class Policy>
  inline OutputIterator bernoulli_series_imp(int start_index,
                                             unsigned number_of_bernoullis_bn,
                                             OutputIterator out_it,
                                             const Policy& pol)
  {
    if(start_index < 0)
    {
       policies::raise_domain_error<T>("boost::math::bernoulli<%1%>", "Start Index should be >= 0 but got %1%", T(start_index), Policy());
    }

    // TBD: Some kind of initialization before main()?
    static std::vector<T> cache_table;

    std::size_t previous_cache_size = cache_table.size();
    std::size_t index = start_index + number_of_bernoullis_bn - 1;

    if(previous_cache_size <= index)
    {
      cache_table.resize(index + 1);

      cache_imp<T,OutputIterator,Policy>(static_cast<int>(previous_cache_size),
                                         static_cast<unsigned int>(index - previous_cache_size + 1),
                                         cache_table.begin() + previous_cache_size,
                                         pol);
    }

    OutputIterator last= out_it + number_of_bernoullis_bn;

    while(out_it!=last)
    {
      *out_it = cache_table[start_index];
      ++out_it;
      ++start_index;
    }

    return out_it;
  }

  template <class T>
  struct max_bernoulli_index
  {
    BOOST_STATIC_CONSTANT(unsigned, value = 17);
  };

  template <>
  struct max_bernoulli_index<float>
  {
    BOOST_STATIC_CONSTANT(unsigned, value = 31);
  };

  template <>
  struct max_bernoulli_index<double>
  {
    BOOST_STATIC_CONSTANT(unsigned, value = 129);
  };

  template <>
  struct max_bernoulli_index<long double>
  {
    BOOST_STATIC_CONSTANT(unsigned, value = 129);
  };

  template <class T>
  inline T unchecked_bernoulli_imp(size_t n, const mpl::int_<3>& )
  {
    static const boost::array<boost::int64_t, 18U> numerators =
    {{
      boost::int64_t(            +1LL),
      boost::int64_t(            +1LL),
      boost::int64_t(            -1LL),
      boost::int64_t(            +1LL),
      boost::int64_t(            -1LL),
      boost::int64_t(            +5LL),
      boost::int64_t(          -691LL),
      boost::int64_t(            +7LL),
      boost::int64_t(         -3617LL),
      boost::int64_t(        +43867LL),
      boost::int64_t(       -174611LL),
      boost::int64_t(       +854513LL),
      boost::int64_t(    -236364091LL),
      boost::int64_t(      +8553103LL),
      boost::int64_t(  -23749461029LL),
      boost::int64_t(+8615841276005LL),
      boost::int64_t(-7709321041217LL),
      boost::int64_t(+2577687858367LL)
    }};

    static const boost::array<boost::int64_t, 18U> denominators =
    {{
      boost::int64_t(      1LL),
      boost::int64_t(      6LL),
      boost::int64_t(     30LL),
      boost::int64_t(     42LL),
      boost::int64_t(     30LL),
      boost::int64_t(     66LL),
      boost::int64_t(   2730LL),
      boost::int64_t(      6LL),
      boost::int64_t(    510LL),
      boost::int64_t(    798LL),
      boost::int64_t(    330LL),
      boost::int64_t(    138LL),
      boost::int64_t(   2730LL),
      boost::int64_t(      6LL),
      boost::int64_t(    870LL),
      boost::int64_t(  14322LL),
      boost::int64_t(    510LL),
      boost::int64_t(      6LL)
    }};
    return T(numerators[n]) / denominators[n];
  }

  template <class T>
  inline T unchecked_bernoulli_imp(size_t n, const mpl::int_<1>& )
  {
    static const boost::array<float, 32U> bernoulli_data =
    {{
      static_cast<T>(+1.00000000000000000000000000000000000000000F),
      static_cast<T>(+0.166666666666666666666666666666666666666667F),
      static_cast<T>(-0.0333333333333333333333333333333333333333333F),
      static_cast<T>(+0.0238095238095238095238095238095238095238095F),
      static_cast<T>(-0.0333333333333333333333333333333333333333333F),
      static_cast<T>(+0.0757575757575757575757575757575757575757576F),
      static_cast<T>(-0.253113553113553113553113553113553113553114F),
      static_cast<T>(+1.16666666666666666666666666666666666666667F),
      static_cast<T>(-7.09215686274509803921568627450980392156863F),
      static_cast<T>(+54.9711779448621553884711779448621553884712F),
      static_cast<T>(-529.124242424242424242424242424242424242424F),
      static_cast<T>(+6192.12318840579710144927536231884057971014F),
      static_cast<T>(-86580.2531135531135531135531135531135531136F),
      static_cast<T>(+1.42551716666666666666666666666666666666667E6F),
      static_cast<T>(-2.72982310678160919540229885057471264367816E7F),
      static_cast<T>(+6.01580873900642368384303868174835916771401E8F),
      static_cast<T>(-1.51163157670921568627450980392156862745098E10F),
      static_cast<T>(+4.29614643061166666666666666666666666666667E11F),
      static_cast<T>(-1.37116552050883327721590879485616327721591E13F),
      static_cast<T>(+4.88332318973593166666666666666666666666667E14F),
      static_cast<T>(-1.92965793419400681486326681448632668144863E16F),
      static_cast<T>(+8.41693047573682615000553709856035437430786E17F),
      static_cast<T>(-4.03380718540594554130768115942028985507246E19F),
      static_cast<T>(+2.11507486380819916056014539007092198581560E21F),
      static_cast<T>(-1.20866265222965259346027311937082525317819E23F),
      static_cast<T>(+7.50086674607696436685572007575757575757576E24F),
      static_cast<T>(-5.03877810148106891413789303052201257861635E26F),
      static_cast<T>(+3.65287764848181233351104308429711779448622E28F),
      static_cast<T>(-2.84987693024508822262691464329106781609195E30F),
      static_cast<T>(+2.38654274996836276446459819192192149717514E32F),
      static_cast<T>(-2.13999492572253336658107447651910973926742E34F),
      static_cast<T>(+2.05009757234780975699217330956723102516667E36F)
    }};

    return bernoulli_data[n];
  }


  template <class T>
  inline T unchecked_bernoulli_imp(size_t n, const mpl::int_<2>& )
  {
    static const boost::array<long double, 130U> bernoulli_data =
    {{
      static_cast<T>(+1.00000000000000000000000000000000000000000L),
      static_cast<T>(+0.166666666666666666666666666666666666666667L),
      static_cast<T>(-0.0333333333333333333333333333333333333333333L),
      static_cast<T>(+0.0238095238095238095238095238095238095238095L),
      static_cast<T>(-0.0333333333333333333333333333333333333333333L),
      static_cast<T>(+0.0757575757575757575757575757575757575757576L),
      static_cast<T>(-0.253113553113553113553113553113553113553114L),
      static_cast<T>(+1.16666666666666666666666666666666666666667L),
      static_cast<T>(-7.09215686274509803921568627450980392156863L),
      static_cast<T>(+54.9711779448621553884711779448621553884712L),
      static_cast<T>(-529.124242424242424242424242424242424242424L),
      static_cast<T>(+6192.12318840579710144927536231884057971014L),
      static_cast<T>(-86580.2531135531135531135531135531135531136L),
      static_cast<T>(+1.42551716666666666666666666666666666666667E6L),
      static_cast<T>(-2.72982310678160919540229885057471264367816E7L),
      static_cast<T>(+6.01580873900642368384303868174835916771401E8L),
      static_cast<T>(-1.51163157670921568627450980392156862745098E10L),
      static_cast<T>(+4.29614643061166666666666666666666666666667E11L),
      static_cast<T>(-1.37116552050883327721590879485616327721591E13L),
      static_cast<T>(+4.88332318973593166666666666666666666666667E14L),
      static_cast<T>(-1.92965793419400681486326681448632668144863E16L),
      static_cast<T>(+8.41693047573682615000553709856035437430786E17L),
      static_cast<T>(-4.03380718540594554130768115942028985507246E19L),
      static_cast<T>(+2.11507486380819916056014539007092198581560E21L),
      static_cast<T>(-1.20866265222965259346027311937082525317819E23L),
      static_cast<T>(+7.50086674607696436685572007575757575757576E24L),
      static_cast<T>(-5.03877810148106891413789303052201257861635E26L),
      static_cast<T>(+3.65287764848181233351104308429711779448622E28L),
      static_cast<T>(-2.84987693024508822262691464329106781609195E30L),
      static_cast<T>(+2.38654274996836276446459819192192149717514E32L),
      static_cast<T>(-2.13999492572253336658107447651910973926742E34L),
      static_cast<T>(+2.05009757234780975699217330956723102516667E36L),
      static_cast<T>(-2.09380059113463784090951852900279701847092E38L),
      static_cast<T>(+2.27526964884635155596492603527692645814700E40L),
      static_cast<T>(-2.62577102862395760473030497361582020814490E42L),
      static_cast<T>(+3.21250821027180325182047923042649852435219E44L),
      static_cast<T>(-4.15982781667947109139170744952623589366896E46L),
      static_cast<T>(+5.69206954820352800238834562191210586444805E48L),
      static_cast<T>(-8.21836294197845756922906534686173330145509E50L),
      static_cast<T>(+1.25029043271669930167323398297028955241772E53L),
      static_cast<T>(-2.00155832332483702749253291988132987687242E55L),
      static_cast<T>(+3.36749829153643742333966769033387530162196E57L),
      static_cast<T>(-5.94709705031354477186604968440515408405791E59L),
      static_cast<T>(+1.10119103236279775595641307904376916046305E62L),
      static_cast<T>(-2.13552595452535011886583850190410656789733E64L),
      static_cast<T>(+4.33288969866411924196166130593792062184514E66L),
      static_cast<T>(-9.18855282416693282262005552155018971389604E68L),
      static_cast<T>(+2.03468967763290744934550279902200200659751E71L),
      static_cast<T>(-4.70038339580357310785752555350060606545967E73L),
      static_cast<T>(+1.13180434454842492706751862577339342678904E76L),
      static_cast<T>(-2.83822495706937069592641563364817647382847E78L),
      static_cast<T>(+7.40642489796788506297508271409209841768797E80L),
      static_cast<T>(-2.00964548027566044834656196727153631868673E83L),
      static_cast<T>(+5.66571700508059414457193460305193569614195E85L),
      static_cast<T>(-1.65845111541362169158237133743199123014950E88L),
      static_cast<T>(+5.03688599504923774192894219151801548124424E90L),
      static_cast<T>(-1.58614682376581863693634015729664387827410E93L),
      static_cast<T>(+5.17567436175456269840732406825071225612408E95L),
      static_cast<T>(-1.74889218402171173396900258776181591451415E98L),
      static_cast<T>(+6.11605199949521852558245252642641677807677E100L),
      static_cast<T>(-2.21227769127078349422883234567129324455732E103L),
      static_cast<T>(+8.27227767987709698542210624599845957312047E105L),
      static_cast<T>(-3.19589251114157095835916343691808148735263E108L),
      static_cast<T>(+1.27500822233877929823100243029266798669572E111L),
      static_cast<T>(-5.25009230867741338994028246245651754469199E113L),
      static_cast<T>(+2.23018178942416252098692981988387281437383E116L),
      static_cast<T>(-9.76845219309552044386335133989802393011669E118L),
      static_cast<T>(+4.40983619784529542722726228748131691918758E121L),
      static_cast<T>(-2.05085708864640888397293377275830154864566E124L),
      static_cast<T>(+9.82144332797912771075729696020975210414919E126L),
      static_cast<T>(-4.84126007982088805087891967099634127611305E129L),
      static_cast<T>(+2.45530888014809826097834674040886903996737E132L),
      static_cast<T>(-1.28069268040847475487825132786017857218118E135L),
      static_cast<T>(+6.86761671046685811921018885984644004360924E137L),
      static_cast<T>(-3.78464685819691046949789954163795568144895E140L),
      static_cast<T>(+2.14261012506652915508713231351482720966602E143L),
      static_cast<T>(-1.24567271371836950070196429616376072194583E146L),
      static_cast<T>(+7.43457875510001525436796683940520613117807E148L),
      static_cast<T>(-4.55357953046417048940633332233212748767721E151L),
      static_cast<T>(+2.86121128168588683453638472510172325229190E154L),
      static_cast<T>(-1.84377235520338697276882026536287854875414E157L),
      static_cast<T>(+1.21811545362210466995013165065995213558174E160L),
      static_cast<T>(-8.24821871853141215484818457296893447301419E162L),
      static_cast<T>(+5.72258779378329433296516498142978615918685E165L),
      static_cast<T>(-4.06685305250591047267679693831158655602196E168L),
      static_cast<T>(+2.95960920646420500628752695815851870426379E171L),
      static_cast<T>(-2.20495225651894575090311752273445984836379E174L),
      static_cast<T>(+1.68125970728895998058311525151360665754464E177L),
      static_cast<T>(-1.31167362135569576486452806355817153004431E180L),
      static_cast<T>(+1.04678940094780380821832853929823089643829E183L),
      static_cast<T>(-8.54328935788337077185982546299082774593270E185L),
      static_cast<T>(+7.12878213224865423522884066771438224721245E188L),
      static_cast<T>(-6.08029314555358993000847118686477458461988E191L),
      static_cast<T>(+5.29967764248499239300942910043247266228490E194L),
      static_cast<T>(-4.71942591687458626443646229013379911103761E197L),
      static_cast<T>(+4.29284137914029810894168296541074669045521E200L),
      static_cast<T>(-3.98767449682322074434477655542938795106651E203L),
      static_cast<T>(+3.78197804193588827138944181161393327898220E206L),
      static_cast<T>(-3.66142336836811912436858082151197348755196E209L),
      static_cast<T>(+3.61760902723728623488554609298914089477541E212L),
      static_cast<T>(-3.64707726451913543621383088655499449048682E215L),
      static_cast<T>(+3.75087554364544090983452410104814189306842E218L),
      static_cast<T>(-3.93458672964390282694891288533713429355657E221L),
      static_cast<T>(+4.20882111481900820046571171111494898242731E224L),
      static_cast<T>(-4.59022962206179186559802940573325591059371E227L),
      static_cast<T>(+5.10317257726295759279198185106496768539760E230L),
      static_cast<T>(-5.78227623036569554015377271242917142512200E233L),
      static_cast<T>(+6.67624821678358810322637794412809363451080E236L),
      static_cast<T>(-7.85353076444504163225916259639312444428230E239L),
      static_cast<T>(+9.41068940670587255245443288258762485293948E242L),
      static_cast<T>(-1.14849338734651839938498599206805592548354E246L),
      static_cast<T>(+1.42729587428487856771416320087122499897180E249L),
      static_cast<T>(-1.80595595869093090142285728117654560926719E252L),
      static_cast<T>(+2.32615353076608052161297985184708876161736E255L),
      static_cast<T>(-3.04957517154995947681942819261542593785327E258L),
      static_cast<T>(+4.06858060764339734424012124124937318633684E261L),
      static_cast<T>(-5.52310313219743616252320044093186392324280E264L),
      static_cast<T>(+7.62772793964343924869949690204961215533859E267L),
      static_cast<T>(-1.07155711196978863132793524001065396932667E271L),
      static_cast<T>(+1.53102008959691884453440916153355334355847E274L),
      static_cast<T>(-2.22448916821798346676602348865048510824835E277L),
      static_cast<T>(+3.28626791906901391668189736436895275365183E280L),
      static_cast<T>(-4.93559289559603449020711938191575963496999E283L),
      static_cast<T>(+7.53495712008325067212266049779283956727824E286L),
      static_cast<T>(-1.16914851545841777278088924731655041783900E290L),
      static_cast<T>(+1.84352614678389394126646201597702232396492E293L),
      static_cast<T>(-2.95368261729680829728014917350525183485207E296L),
      static_cast<T>(+4.80793212775015697668878704043264072227967E299L),
      static_cast<T>(-7.95021250458852528538243631671158693036798E302L),
      static_cast<T>(+1.33527841873546338750122832017820518292039E306L)
    }};

    return bernoulli_data[n];
  }

  template<class T>
  inline T unchecked_bernoulli_b2n(size_t n)
  {
    typedef typename mpl::if_c<
         (std::numeric_limits<T>::max_exponent >= 128)
      && (std::numeric_limits<T>::max_exponent < 1024)
      && (std::numeric_limits<T>::digits10 <= std::numeric_limits<float>::digits10)
      && boost::is_convertible<T, float>::value,

      mpl::int_<1>,

      typename mpl::if_c<
           (std::numeric_limits<T>::max_exponent >= 1024)
        && (std::numeric_limits<T>::digits10 <= std::numeric_limits<long double>::digits10)
        && boost::is_convertible<T, long double>::value,
        mpl::int_<2>,
        mpl::int_<3>
      >::type
    >::type tag_type;

    return unchecked_bernoulli_imp<T>(n, tag_type());
  }

} } } // namespace boost::math::detail

#endif // _BOOST_BERNOULLI_B2N_2013_05_30_HPP_
