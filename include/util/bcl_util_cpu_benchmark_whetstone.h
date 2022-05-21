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

#ifndef BCL_UTIL_CPU_BENCHMARK_WHETSTONE_H_
#define BCL_UTIL_CPU_BENCHMARK_WHETSTONE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"
#include "bcl_util_stopwatch.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CPUBenchmarkWhetstone
    //! @brief this class performs a cpu benchmark according to the established whetstone procedure
    //! @details the benchmark is performed for a given duration, and the result are operations or floating point
    //!          operations per second, that are stored in a table for various types of operations
    //!
    //! @tparam t_DataType the data type float, double for which the benchmark is performed
    //!
    //! C/C++ Whetstone Benchmark Single or Double Precision
    //!
    //! Original concept        Brian Wichmann NPL      1960's
    //! Original author         Harold Curnow  CCTA     1972
    //! Self timing versions    Roy Longbottom CCTA     1978/87
    //! Optimisation control    Bangor University       1987/90
    //! C/C++ Version           Roy Longbottom          1996
    //! Compatibility & timers  Al Aburto               1996
    //!
    //!          Official version approved by:
    //!     Harold Curnow  100421.1615@compuserve.com
    //!  Happy 25th birthday Whetstone, 21 November 1997
    //!
    //! @see http://www.roylongbottom.org.uk/whetstone.pdf
    //! @see @link example_util_cpu_benchmark_whetstone.cpp @endlink
    //! @author woetzen
    //! @date Oct 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class CPUBenchmarkWhetstone :
      public ObjectInterface
    {

      enum SectionType
      {
        e_N1FloatingPoint,     //!< FloatingPoint1
        e_N2FloatingPoint,     //!< FloatingPoint2
        e_N3Condition,         //!< Condition
        e_N4Integer,           //!< Integer
        e_N5Trigonometry,      //!< Trigonometry
        e_N6FloatingPoint,     //!< FloatingPoint3
        e_N7Assignments,       //!< Assignments
        e_N8ExponentialSquare, //!< ExponentialSquare
        e_Sum,
        e_Average,
        s_NumberSectionTypes
      };

      enum ResultType
      {
        e_ResultOperations,              //!< Result
        e_ResultTime,                    //!< Seconds
        e_ResultMegaOperationsPerSecond, //!< mega = 10^6 operations per second
        s_NumberResultTypes
      };

    private:

    //////////
    // data //
    //////////

      Time m_TimeUsed; //!< times used for the benchmark
      t_DataType m_Check;    //!< will be used to check if test did run
      storage::Table< t_DataType> m_ResultTable;

      static const size_t *GetSectionMultiplier();
      static const std::string *GetLoopTableHeader();
      static const std::string *GetRowNames();

    public:

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CPUBenchmarkWhetstone() :
        m_Check( 0.0),
        m_ResultTable( storage::TableHeader( storage::Vector< std::string>( s_NumberResultTypes, GetLoopTableHeader())))
      {
      }

      //! @brief constructor from table
      //! @param RESULT_TABLE result table
      CPUBenchmarkWhetstone( const storage::Table< t_DataType> &RESULT_TABLE) :
        m_Check( 0.0),
        m_ResultTable( RESULT_TABLE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new CPUBenchmarkWhetstone
      CPUBenchmarkWhetstone *Clone() const
      {
        return new CPUBenchmarkWhetstone< t_DataType>();
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief access result table
      //! @return reference to results
      const storage::Table< t_DataType> &GetResultTable() const
      {
        return m_ResultTable;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief run the benchmark
      //! @param DURATION the time the benchmark should run, the longer the more accurate it will be, 100 seconds is recommended
      //! @return 0 on success, -1 if the benchmark did not work
      int RunBenchmark( const Time &DURATION)
      {
        // empty the result table
        m_ResultTable.Reset();

        int count = 10, calibrate = 1;
        size_t xtra = 1;
        const size_t x100( 100);

        BCL_MessageStd( "Calibrate");
        do
        {
          m_TimeUsed.SetZero();

          RunWhetstonesTests( xtra, x100, calibrate);

          BCL_MessageStd( Format()( m_TimeUsed.GetSecondsFractional()) + " Seconds " + util::Format()( ( t_DataType)( xtra)) + " Passes (x 100)");
          calibrate = calibrate + 1;
          count = count - 1;
          if( m_TimeUsed > Time( 2, 0))
            count = 0;
          else
            xtra = xtra * 5;
        }
        while( count > 0);
        if( m_TimeUsed > Time())
          xtra = size_t( ( DURATION.GetSecondsFractional() * xtra) / m_TimeUsed.GetSecondsFractional());
        if( xtra < 1) xtra = 1;

        calibrate = 0;
        BCL_MessageStd( "Use " + util::Format()( xtra) + "  passes (x 100)");

        m_TimeUsed.SetZero();
        RunWhetstonesTests( xtra, x100, calibrate);

        t_DataType mwips;
        if( m_TimeUsed > Time())
          mwips = t_DataType( xtra) * t_DataType( x100) / t_DataType( 10 * m_TimeUsed.GetSecondsFractional());
        else
          mwips = 0;

        BCL_MessageStd( "MWIPS            " + util::Format()( mwips) + " " + util::Format()( m_TimeUsed.GetSecondsFractional()));
        if( m_Check == 0)
        {
          BCL_MessageStd( "Wrong answer");
          m_ResultTable.Reset();
          return -1;
        }

        AddSumToTable();

        return 0;
      }

      //! @brief get average operations per second
      //! @brief calculate the average operations per second
      t_DataType GetAverageOperationsPerSecond() const
      {
        return m_ResultTable[ GetRowNames()[ e_Average]][ GetLoopTableHeader()[ e_ResultMegaOperationsPerSecond]] * t_DataType( 1000000);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief run all whetstones tests
      //! @param XTRA the factor for additional iterations, which is used to adjust the test for the target duration time
      //! @param X100 factor of 100 which is the basic multiplier of the number of iterations
      //! @param CALIBRATE if non zero, it is a calibration run to aim for the target duration and no result is generated
      void RunWhetstonesTests( const size_t XTRA, const size_t X100, int CALIBRATE)
      {
        // stopwatch for timing
        Stopwatch stop_watch( false);

        // variables
        size_t i, ix, n_max;
        t_DataType x, y, z;
        long j, k, l;
        t_DataType e1[ 4];

        t_DataType t  = t_DataType( 0.49999975);
        t_DataType t0 = t;
        const t_DataType t1( t_DataType( 0.50000025));
        const t_DataType t2( t_DataType( 2.0));

        // set the checksum to 0
        m_Check = 0.0;

        // Section 1, Array elements
        e1[ 0] =  1.0;
        e1[ 1] = -1.0;
        e1[ 2] = -1.0;
        e1[ 3] = -1.0;
        n_max = GetSectionMultiplier()[ e_N1FloatingPoint] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            for( i = 0; i < n_max; ++i)
            {
              e1[ 0] = ( e1[ 0] + e1[ 1] + e1[ 2] - e1[ 3]) * t;
              e1[ 1] = ( e1[ 0] + e1[ 1] - e1[ 2] + e1[ 3]) * t;
              e1[ 2] = ( e1[ 0] - e1[ 1] + e1[ 2] + e1[ 3]) * t;
              e1[ 3] = (-e1[ 0] + e1[ 1] + e1[ 2] + e1[ 3]) * t;
            }
            t = t_DataType( 1.0) - t;
          }
          t = t0;
        }
        stop_watch.Stop();
        AddToResultTable( e_N1FloatingPoint, t_DataType( n_max) * t_DataType( 16.0) * t_DataType( XTRA), e1[ 3], stop_watch.GetTotalTime(), CALIBRATE);

        // Section 2, Array as parameter
        n_max = GetSectionMultiplier()[ e_N2FloatingPoint] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            for( i = 0; i < n_max; ++i)
            {
              PA( e1, t, t2);
            }
            t = t_DataType( 1.0) - t;
          }
          t = t0;
        }
        stop_watch.Stop();
        AddToResultTable( e_N2FloatingPoint, t_DataType( n_max) * t_DataType( 96.0) * t_DataType( XTRA), e1[ 3], stop_watch.GetTotalTime(), CALIBRATE);

        // Section 3, Conditional jumps
        j = 1;
        n_max = GetSectionMultiplier()[ e_N3Condition] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            for( i = 0; i < n_max; ++i)
            {
              if( j == 1) j = 2;
              else        j = 3;
              if( j > 2)  j = 0;
              else        j = 1;
              if( j < 1)  j = 1;
              else        j = 0;
            }
          }
        }
        stop_watch.Stop();
        AddToResultTable( e_N3Condition, t_DataType( n_max) * t_DataType( 3.0) * t_DataType( XTRA), t_DataType( j), stop_watch.GetTotalTime(), CALIBRATE);

        // Section 4, Integer arithmetic
        j = 1;
        k = 2;
        l = 3;
        n_max = GetSectionMultiplier()[ e_N4Integer] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            for( i = 0; i < n_max; ++i)
            {
              j = j * ( k - j) * ( l - k);
              k = l * k - ( l - j) * k;
              l = ( l - k) * ( k + j);
              e1[ l - 2] = t_DataType( j + k + l);
              e1[ k - 2] = t_DataType( j * k * l);
            }
          }
        }
        stop_watch.Stop();
        x = e1[ 0] + e1[ 1];
        AddToResultTable( e_N4Integer, t_DataType( n_max) * t_DataType( 15.0) * t_DataType( XTRA), x, stop_watch.GetTotalTime(), CALIBRATE);

        // Section 5, Trig functions
        x = 0.5;
        y = 0.5;
        n_max = GetSectionMultiplier()[ e_N5Trigonometry] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            for( i = 1; i < n_max; ++i)
            {
              x = t * atan( t2 * sin( x) * cos( x) / ( cos( x + y) + cos( x - y) - t_DataType( 1.0)));
              y = t * atan( t2 * sin( y) * cos( y) / ( cos( x + y) + cos( x - y) - t_DataType( 1.0)));
            }
            t = t_DataType( 1.0) - t;
          }
          t = t0;
        }
        stop_watch.Stop();
        AddToResultTable( e_N5Trigonometry, t_DataType( n_max) * t_DataType( 26.0) * t_DataType( XTRA), y, stop_watch.GetTotalTime(), CALIBRATE);

        // Section 6, Procedure calls
        x = 1.0;
        y = 1.0;
        z = 1.0;
        n_max = GetSectionMultiplier()[ e_N6FloatingPoint] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            for( i = 0; i < n_max; ++i)
            {
              P3( &x, &y, &z, t, t1, t2);
            }
          }
        }
        stop_watch.Stop();
        AddToResultTable( e_N6FloatingPoint, t_DataType( n_max) * t_DataType( 6.0) * t_DataType( XTRA), z, stop_watch.GetTotalTime(), CALIBRATE);

        // Section 7, Array refrences
        j = 0;
        k = 1;
        l = 2;
        e1[ 0] = 1.0;
        e1[ 1] = 2.0;
        e1[ 2] = 3.0;
        n_max = GetSectionMultiplier()[ e_N7Assignments] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            for( i = 0; i < n_max; ++i)
            {
              PO( e1, j, k, l);
            }
          }
        }
        stop_watch.Stop();
        AddToResultTable( e_N7Assignments, t_DataType( n_max) * t_DataType( 3.0) * t_DataType( XTRA), e1[ 2], stop_watch.GetTotalTime(), CALIBRATE);

        // Section 8, Standard functions
        n_max = GetSectionMultiplier()[ e_N8ExponentialSquare] * X100;
        stop_watch.Reset();
        stop_watch.Start();
        {
          for( ix = 0; ix < XTRA; ++ix)
          {
            x = 0.75; // if XTRA is too large, x becomes one, which seems to be absolutely slow on any cpu
            for( i = 0; i < n_max; ++i)
            {
              x = sqrt( exp( log( x) / t1));
            }
          }
        }
        stop_watch.Stop();
        AddToResultTable( e_N8ExponentialSquare, t_DataType( n_max) * t_DataType( 4.0) * t_DataType( XTRA), x, stop_watch.GetTotalTime(), CALIBRATE);
      }

      //! @brief add result to table
      //! @param SECTION_TYPE which section of the whetstone test was performed
      //! @param NUMBER_OPERATIONS the number of operations that were performed
      //! @param checknum the check num to verify the success fo the test
      //! @param TIME the time this section ran for
      //! @param CALIBRATE if non zero, it was a calibration run, else result is added to the result table
      void AddToResultTable( const SectionType SECTION_TYPE, t_DataType NUMBER_OPERATIONS, t_DataType checknum, const Time &TIME, int CALIBRATE)
      {
        m_Check    += checknum;
        m_TimeUsed += TIME;

        // do not insert results for calibration rounds
        if( CALIBRATE != 0)
        {
          return;
        }

        t_DataType operations_per_second( 0);

        // if time has passed
        if( TIME > Time())
        {
          operations_per_second = NUMBER_OPERATIONS / t_DataType( 1000000L * TIME.GetSecondsFractional());
        }

        // insert into result table
        storage::Row< t_DataType> &row( m_ResultTable.InsertRow( GetRowNames()[ SECTION_TYPE]));
        row( e_ResultOperations)              = NUMBER_OPERATIONS;
        row( e_ResultTime)                    = t_DataType( TIME.GetSecondsFractional());
        row( e_ResultMegaOperationsPerSecond) = operations_per_second;
      }

      //! @brief sums up all iterations and adds that to the table
      void AddSumToTable()
      {
        // iterate over all rows, add up total time and operations
        t_DataType operations( 0);
        t_DataType time( 0);

        for( typename storage::Table< t_DataType>::const_iterator itr( m_ResultTable.Begin()), itr_end( m_ResultTable.End()); itr != itr_end; ++itr)
        {
          operations += itr->Second()( e_ResultOperations);
          time       += itr->Second()( e_ResultTime);
        }

        // insert sum row
        storage::Row< t_DataType> &row_sum( m_ResultTable.InsertRow( GetRowNames()[ e_Sum]));
        row_sum( e_ResultOperations             ) = operations;
        row_sum( e_ResultMegaOperationsPerSecond) = GetUndefined< t_DataType>();
        row_sum( e_ResultTime                   ) = time;

        // insert average row
        storage::Row< t_DataType> &row_average( m_ResultTable.InsertRow( GetRowNames()[ e_Average]));
        row_average( e_ResultOperations             ) = GetUndefined< t_DataType>();
        row_average( e_ResultMegaOperationsPerSecond) = operations / t_DataType( 1000000) / time;
        row_average( e_ResultTime                   ) = GetUndefined< t_DataType>();
      }

      //! @brief helper function for floating point operations 2
      void PA( t_DataType e[ 4], t_DataType t, t_DataType t2)
      {
        for( size_t j( 0); j < 6; ++j)
        {
          e[ 0] = ( e[ 0] + e[ 1] + e[ 2] - e[ 3]) * t;
          e[ 1] = ( e[ 0] + e[ 1] - e[ 2] + e[ 3]) * t;
          e[ 2] = ( e[ 0] - e[ 1] + e[ 2] + e[ 3]) * t;
          e[ 3] = (-e[ 0] + e[ 1] + e[ 2] + e[ 3]) / t2;
        }
      }

      //! helper function for assignment operations
      void PO( t_DataType e1[ 4], long j, long k, long l)
      {
        e1[ j] = e1[ k];
        e1[ k] = e1[ l];
        e1[ l] = e1[ j];
      }

      //! @brief helper function for floating point operations 3
      void P3( t_DataType *x, t_DataType *y, t_DataType *z, t_DataType t, t_DataType t1, t_DataType t2)
      {
        *x = *y;
        *y = *z;
        *x = t  * ( *x + *y);
        *y = t1 * ( *x + *y);
        *z = ( *x + *y) / t2;
      }

    }; // template class CPUBenchmarkWhetstone

    //! @brief instance of s_SectionMultiplier
    template< typename t_DataType>
    const size_t *CPUBenchmarkWhetstone< t_DataType>::GetSectionMultiplier()
    {
      static const size_t s_section_multiplier[] =
      {
        120,  // e_N1FloatingPoint
        14,   // e_N2FloatingPoint
        345,  // e_N3Condition
        210,  // e_N4Integer
        32,   // e_N5Trigonometry
        899,  // e_N6FloatingPoint
        616,  // e_N7Assignments
        93    // e_N8ExponentialSquare
      };

      // end
      return s_section_multiplier;
    }

    //! @brief instance of s_LoopTableHeader
    template< typename t_DataType>
    const std::string *CPUBenchmarkWhetstone< t_DataType>::GetLoopTableHeader()
    {
      static const std::string s_loop_table_header[] =
      {
        "TotalOperations",
        "Time[s]",
        "OperationsPerSecond"
      };

      // end
      return s_loop_table_header;
    }

    //! @brief instance of s_RowNames
    template< typename t_DataType>
    const std::string *CPUBenchmarkWhetstone< t_DataType>::GetRowNames()
    {
      static const std::string s_row_names[] =
      {
        "FloatingPoint1",    // e_N1FloatingPoint
        "FloatingPoint2",    // e_N2FloatingPoint
        "Condition",         // e_N3Condition
        "Integer",           // e_N4Integer
        "Trigonometry",      // e_N5Trigonometry
        "FloatingPoint3",    // e_N6FloatingPoint
        "Assignments",       // e_N7Assignments
        "ExponentialSquare", // e_N8ExponentialSquare
        "Sum",               // e_Sum
        "Average"            // e_Average
      };

      // end
      return s_row_names;
    }

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_CPU_BENCHMARK_WHETSTONE_H_
