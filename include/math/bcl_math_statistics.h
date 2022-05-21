// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_MATH_STATISTICS_H_
#define BCL_MATH_STATISTICS_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "random/bcl_random_uniform_distribution.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

// disable warning c4996 for use of unchecked iterators in stl algortihm functions
#if defined (_MSC_VER)
  #pragma warning ( disable : 4996)
#endif

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Statistics
    //! @brief collects methods to determine statistical values like min/max, mean, norm on iterator ranges
    //!
    //! @see @link example_math_statistics.cpp @endlink
    //! @author mueller
    //! @date Feb 11, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Statistics
    {
    private:

      Statistics()
      {
      }

    public:

      // these functions iterate over a range to return a result
      // For convenience, these functions all use typename std::iterator_traits< t_InputIterator>::value_type,
      // which gives the type that the iterator returns when dereferenced.  This way, the user does not have to
      // specify the value type of iterator as a template parameter

      //! @brief get the minimum value of a sequence
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        MinimumValue( t_InputIterator BEGIN, t_InputIterator END)
      {
        return *std::min_element( BEGIN, END);
      }

      //! @brief get a pointer to the minimum value of a sequence
      template< typename t_InputIterator>
      static t_InputIterator PointerToMinimumValue( t_InputIterator BEGIN, t_InputIterator END)
      {
        return std::min_element( BEGIN, END);
      }

      //! @brief get the index of the minimum value of a sequence
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::difference_type
        MinimumIndex( t_InputIterator BEGIN, t_InputIterator END)
      {
        return std::distance( BEGIN, std::min_element( BEGIN, END));
      }

      //! @brief get the maximum value of a sequence
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        MaximumValue( t_InputIterator BEGIN, t_InputIterator END)
      {
        return *std::max_element( BEGIN, END);
      }

      //! @brief get a pointer to the maximum value of a sequence
      template< typename t_InputIterator>
      static t_InputIterator PointerToMaximumValue( t_InputIterator BEGIN, t_InputIterator END)
      {
        return std::max_element( BEGIN, END);
      }

      //! @brief get the index of the maximum value
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::difference_type
        MaximumIndex( t_InputIterator BEGIN, t_InputIterator END)
      {
        return std::distance( BEGIN, std::max_element( BEGIN, END));
      }

      //! @brief get the sum of a sequence
      //! @param BEGIN, END iterators on a sequence
      //! @param START_VALUE initial value for the sum, defaults to 0
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        Sum( const t_InputIterator BEGIN, const t_InputIterator END,
             const typename std::iterator_traits< t_InputIterator>::value_type &START_VALUE
               = ( typename std::iterator_traits< t_InputIterator>::value_type)( 0))
      {
        return std::accumulate( BEGIN, END, START_VALUE);
      }

      //! compute and return square of norm
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        SquareNorm( t_InputIterator BEGIN, const t_InputIterator END)
      {
        return std::inner_product
               (
                 BEGIN,
                 END,
                 BEGIN,
                 ( typename std::iterator_traits< t_InputIterator>::value_type)( 0)
               );
      }

      //! compute and return Norm of VectorMatrixTensorBase
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        Norm( t_InputIterator BEGIN, const t_InputIterator END)
      {
        return Sqrt( SquareNorm( BEGIN, END));
      }

      //! compute the mean of all elements
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        Mean( const t_InputIterator BEGIN, const t_InputIterator END)
      {
        typename std::iterator_traits< t_InputIterator>::value_type size( END - BEGIN);
        return Sum( BEGIN, END) / size;
      }

      //! compute the variance of all elements
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        Variance( t_InputIterator BEGIN, const t_InputIterator END)
      {
        typename std::iterator_traits< t_InputIterator>::value_type size( END - BEGIN);
        return std::abs( SquareNorm( BEGIN, END) / size - Sqr( Mean( BEGIN, END)));
      }

      //! compute the variance of all elements
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type
        StandardDeviation( t_InputIterator BEGIN, const t_InputIterator END)
      {
        return Sqrt( Variance( BEGIN, END));
      }

      //! normalize all elements
      template< typename t_InputIterator>
      static void Normalize( t_InputIterator BEGIN, const t_InputIterator END)
      {
        // calculate norm
        const typename std::iterator_traits< t_InputIterator>::value_type norm( Norm( BEGIN, END));

        //if the norm is zero then the vector will not be divided by zero
        if( norm == typename std::iterator_traits< t_InputIterator>::value_type( 0))
        {
          BCL_MessageDbg( "returning since norm is " + util::Format()( norm));
          return;
        }

        // divide each element in the range by the norm
        std::transform
        (
          BEGIN,
          END,
          BEGIN,
          std::bind2nd( std::divides< typename std::iterator_traits< t_InputIterator>::value_type>(), norm)
        );
      }

      //! @brief calculate the correlation coefficient (r) by taking the square-root of the coefficient of determination (r^2)
      //! @param BEGIN_A iterator of container A pointing to beginning
      //! @param END_A iterator of container A pointing to end
      //! @param BEGIN_B iterator of container B pointing to beginning
      //! @param END_B iterator of container B pointing to end
      //! @return pearson correlation coeffcient
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type CorrelationCoefficient
      (
        t_InputIterator BEGIN_A, const t_InputIterator END_A,
        t_InputIterator BEGIN_B, const t_InputIterator END_B
      )
      {
        return Sqrt( RSquared( BEGIN_A, END_A, BEGIN_B, END_B));
      }

      //! calculate the correlation as cor( x, y) = covariance( x, y) / (stddev( x) * stddev( y)) (Pearson correlation coefficient)
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type Correlation
      (
        t_InputIterator BEGIN_A, const t_InputIterator END_A,
        t_InputIterator BEGIN_B, const t_InputIterator END_B
      )
      {
        return
            Covariance( BEGIN_A, END_A, BEGIN_B, END_B)
            / ( StandardDeviation( BEGIN_A, END_A) * StandardDeviation( BEGIN_B, END_B));
      }

      //! @brief calculate the Pearson correlation coefficient between the ranked variables (Spearman correlation coefficient)
      //! @param BEGIN_A iterator of container A pointing to beginning
      //! @param END_A iterator of container A pointing to end
      //! @param BEGIN_B iterator of container B pointing to beginning
      //! @param END_B iterator of container B pointing to end
      //! @return spearman correlation coeffcient
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type CorrelationSpearman
      (
        const t_InputIterator BEGIN_A, const t_InputIterator END_A,
        const t_InputIterator BEGIN_B, const t_InputIterator END_B
      )
      {
        // check input containers for same size
        BCL_Assert( std::distance( BEGIN_A, END_A) == std::distance( BEGIN_B, END_B), "CorrelationSpearman: inputs have to have same size!");

        // maps containing the values of containers as key and occurences as value
        // the automatic sorting of a map takes care of the ranking order of elements
        storage::Map< typename std::iterator_traits< t_InputIterator>::value_type, double> ranks_map_a, ranks_map_b;

        // initialize iterators
        t_InputIterator itr_a( BEGIN_A);
        t_InputIterator itr_b( BEGIN_B);

        // add values of container and update their respective count
        for( ; itr_a != END_A; ++itr_a, ++itr_b)
        {
          ranks_map_a[ *itr_a] += 1.0;
          ranks_map_b[ *itr_b] += 1.0;
        }

        // counter for ranks of container a
        double rank_a( 0);

        // iterate over all values in rank map a and determine the respective rank
        // elements with same rank share the average rank
        for
        (
          typename storage::Map< typename std::iterator_traits< t_InputIterator>::value_type, double>::iterator
            itr_map_a( ranks_map_a.Begin()), itr_map_a_end( ranks_map_a.End());
          itr_map_a != itr_map_a_end;
          ++itr_map_a
        )
        {
          // count for key in rank map a
          double count( itr_map_a->second);
          // adjust rank in map and compute average rank
          itr_map_a->second = rank_a + ( count + 1) / 2.0;
          // increment rank
          rank_a += count;
        }

        // counter for ranks of container b
        double rank_b( 0);

        // iterate over all values in rank map b and determine the respective rank
        // elements with same rank share the average rank
        for
        (
          typename storage::Map< typename std::iterator_traits< t_InputIterator>::value_type, double>::iterator
            itr_map_b( ranks_map_b.Begin()), itr_map_b_end( ranks_map_b.End());
          itr_map_b != itr_map_b_end;
          ++itr_map_b
        )
        {
          // count for key in rank map b
          double count( itr_map_b->second);
          // adjust rank in map and compute average rank
          itr_map_b->second = rank_b + ( count + 1) / 2.0;
          // increment rank
          rank_b += count;
        }

        // initialize container for rankings
        storage::Vector< double> ranks_a, ranks_b;

        // fill final ranking containers
        for( itr_a = BEGIN_A, itr_b = BEGIN_B; itr_a != END_A; ++itr_a, ++itr_b)
        {
          ranks_a.PushBack( ranks_map_a[ *itr_a]);
          ranks_b.PushBack( ranks_map_b[ *itr_b]);
        }

        // calculate the Pearson correlation coefficient between the ranked variables (Spearman)
        return Correlation( ranks_a.Begin(), ranks_a.End(), ranks_b.Begin(), ranks_b.End());
      }

      //!calculate the coefficient of determination between x and y (r^2) with the following formula
      //!R^2 = cov(x,y)^2/(var(x)*var(y))-
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type RSquared
      (
        t_InputIterator BEGIN_A, const t_InputIterator END_A,
        t_InputIterator BEGIN_B, const t_InputIterator END_B
      )
      {
        double covariance( Covariance( BEGIN_A, END_A, BEGIN_B, END_B));
        double variance_a( Variance( BEGIN_A, END_A));
        double variance_b( Variance( BEGIN_B, END_B));
        return Sqr( covariance) / ( variance_a * variance_b);
      }

      //!calculate the covariance between x and y with the following formula
      //!Cov(x,y) (mean of the products x*y) - (mean of x) * (mean of y)
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type Covariance
      (
        t_InputIterator BEGIN_A, const t_InputIterator END_A,
        t_InputIterator BEGIN_B, const t_InputIterator END_B
      )
      {
        double products( ScalarProduct( BEGIN_A, END_A, BEGIN_B, END_B));
        return products / ( END_A - BEGIN_A) - Mean( BEGIN_A, END_A) * Mean( BEGIN_B, END_B);
      }

      //!calculate the covariance between x and y with the following formula
      //!Cov(x,y) (mean of the products x*y) - (mean of x) * (mean of y)
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type ScalarProduct
      (
        t_InputIterator BEGIN_A, const t_InputIterator END_A,
        t_InputIterator BEGIN_B, const t_InputIterator END_B
      )
      {
        BCL_Assert
        (
          std::distance( BEGIN_A, END_A) == std::distance( BEGIN_B, END_B),
          "cannot calculate the ScalarProduct between data sets of different sizes"
        );

        return
          std::inner_product
          (
            BEGIN_A,
            END_A,
            BEGIN_B,
            ( typename std::iterator_traits< t_InputIterator>::value_type)( 0)
          );
      }

      //! set all elements in array to random numbers
      template< typename t_DataType, typename t_InputIterator>
      static void SetRand
      (
        t_InputIterator BEGIN,
        const t_InputIterator END,
        const t_DataType MIN,
        const t_DataType MAX
      )
      {
        //randomize each element
        for( ; BEGIN != END; ++BEGIN)
        {
          ( *BEGIN) = random::GetGlobalRandom().Random< t_DataType>( MIN, MAX);
        }
      }

      //! set to sum SUM (sum of all elements = SUM)
      template< typename t_DataType, typename t_InputIterator>
      static void SetToSum
      (
        t_InputIterator BEGIN,
        const t_InputIterator END,
        const t_DataType SUM = 1
      )
      {
        BCL_Assert( SUM != t_DataType( 0), "passing t_DataType( 0) does not make sense!");
        const t_DataType this_sum( Sum( BEGIN, END));
        if( this_sum == t_DataType( 0))
        {
          return;
        }

        // multiply each with the fraction of the target SUM vs. the current sum of that range
        std::transform( BEGIN, END, BEGIN, std::bind2nd( std::multiplies< t_DataType>(), SUM / this_sum));
      }

      //! checks that every position is not of type undefined
      template< typename t_InputIterator>
      static bool IsDefined( t_InputIterator BEGIN, const t_InputIterator END)
      {
        //iterate over all elements
        for( ; BEGIN != END; ++BEGIN)
        {
          //is current pointer value of undefined type (double - nan, int - util::GetUndefined< int>())
          if( !util::IsDefined( *BEGIN))
          {
            return false;
          }
        }

        //end
        return true;
      }

      //! Cumulative Euclidian
      //! sums the difference between the sum of increasingly large portions of the ranges
      template< typename t_InputIterator>
      static typename std::iterator_traits< t_InputIterator>::value_type CumulativeEuclidian
      (
        t_InputIterator BEGIN_A, const t_InputIterator END_A,
        t_InputIterator BEGIN_B, const t_InputIterator END_B
      )
      {
        // create values to hold the total and cumulative of each range
        typename std::iterator_traits< t_InputIterator>::value_type
          cumulative_a( 0), cumulative_b( 0), cumulative_distance_squared( 0);

        // iterate over the two iterator ranges
        for
        (
          t_InputIterator itr_a( BEGIN_A), itr_b( BEGIN_B);
          itr_a != END_A && itr_b != END_B;
          ++itr_a, ++itr_b
        )
        {
          cumulative_a += *itr_a;
          cumulative_b += *itr_b;
          cumulative_distance_squared += Sqr( cumulative_a - cumulative_b);
        }
        return Sqrt( cumulative_distance_squared);
      }

      //! @brief calculates the critical value for a given KS D value
      //! @param ALPHA the confidence level
      //! @param N1 one sample size
      //! @param N2 the other sample size
      //! @return the critical value; this is a table lookup, so will return the critical value for the next-highest alpha value
      //!         if the actual alpha value is not in the table (e.g. uncommon alpha value)
      static double GetKolmogorovSmirnovCriticalValue( const double &ALPHA, const size_t &N1, const size_t &N2)
      {
        static double alpha_table[] = { 0.1, 0.05, 0.025, 0.01, 0.005, 0.001};
        static double crit_table[] = { 1.22, 1.36, 1.48, 1.63, 1.73, 1.95};

        if( N1 == 0 || N2 == 0)
        {
          return util::GetUndefined< double>();
        }

        size_t index( 0);
        for
        ( 
          ; index < 6 
            && !EqualWithinAbsoluteTolerance( ALPHA, alpha_table[ index], 0.0001) 
            && ALPHA < alpha_table[ index]; 
          ++index
        );
        double coeff( crit_table[ index]);
        return coeff * std::sqrt( double( N1 + N2) / double( N1 * N2));
      }

      //! @brief Kolmogorov-Smirnov test; determines if two samples are from the same distribution
      //! @return the KS D and critical value (see https://en.wikipedia.org/wiki/Kolmogorov¿Smirnov_test)
      template< typename t_InputIterator>
      static std::pair< double, double> KolmogorovSmirnov
      (
        const t_InputIterator &BEGIN_A, const t_InputIterator &END_A,
        const t_InputIterator &BEGIN_B, const t_InputIterator &END_B,
        const double &ALPHA = 0.05
      )
      {
        typedef typename std::iterator_traits< t_InputIterator>::value_type val_type;

        size_t a_size( std::distance( BEGIN_A, END_A)), b_size( std::distance( BEGIN_B, END_B));
        if( a_size == 0 || b_size == 0)
        {
          return std::pair< double, double>( util::GetUndefined< double>(), util::GetUndefined< double>());
        }
        
        // generate vectors from the iterators, then sort them in descending order
        std::vector< val_type> dist_a( BEGIN_A, END_A);
        std::vector< val_type> dist_b( BEGIN_B, END_B);
        std::sort( dist_a.begin(), dist_a.end(), std::less< val_type>());
        std::sort( dist_b.begin(), dist_b.end(), std::less< val_type>());

        // current KS statistic
        double ks( 0);

        typename std::vector< val_type>::const_iterator itr_a( dist_a.begin()), itr_a_end( dist_a.end());
        typename std::vector< val_type>::const_iterator itr_b( dist_b.begin()), itr_b_end( dist_b.end());

        size_t n_a( 0), n_b( 0);
        val_type val( std::min( *itr_a, *itr_b));
        while( itr_a != itr_a_end && itr_b != itr_b_end)
        {
          // move itr_a to the first position greater than val
          for( ; itr_a != itr_a_end && *itr_a <= val; ++n_a, ++itr_a);

          // move itr_b to the first position greater than val
          for( ; itr_b != itr_b_end && *itr_b <= val; ++n_b, ++itr_b);

          double p_a( double( n_a) / a_size);
          double p_b( double( n_b) / b_size);
          double diff( std::fabs( p_a - p_b));
          if( diff > ks)
          {
            ks = diff;
          }

          if( itr_a != itr_a_end && itr_b != itr_b_end)
          {
            val = std::min( *itr_a, *itr_b);
          }
          else if( itr_b == itr_b_end)
          {
            val = *itr_a;
          }
          else
          {
            val = *itr_b;
          }
        }

        return std::pair< double, double>( ks, GetKolmogorovSmirnovCriticalValue( ALPHA, a_size, b_size));
      }

    }; // class Statistics

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_STATISTICS_H_
