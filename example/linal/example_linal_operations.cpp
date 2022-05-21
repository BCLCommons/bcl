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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "linal/bcl_linal_operations.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_si_ptr_vector.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_operations.cpp
  //!
  //! @author woetzen, loweew
  //! @date Mar 15, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalOperations :
    public ExampleInterface
  {
    // class generator:
    template< typename t_DataType>
    struct Increment
    {
      t_DataType m_Current;
      Increment() :
        m_Current( 0)
      {
      }

      t_DataType operator()()
      {
        return m_Current += t_DataType( 0.1);
      }

    }; // struct Increment

  public:

    ExampleLinalOperations *Clone() const
    {
      return new ExampleLinalOperations( *this);
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

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if all vectors are equal
    //! @param RESULTS map of operations string and the result vector to be compare to other vectors
    //! @return true, if all vectors are equal within a small margin
    template< typename t_DataType>
    bool EqualWithinTolerance( const storage::Map< std::string, linal::Vector< t_DataType> > &RESULTS) const
    {
      bool equal( true);

      t_DataType deviation( 0);
      // iterate over all result pairs
      for( typename storage::Map< std::string, linal::Vector< t_DataType> >::const_iterator itr1( RESULTS.Begin()), itr_end( RESULTS.End()); itr1 != itr_end; ++itr1)
      {
        typename storage::Map< std::string, linal::Vector< t_DataType> >::const_iterator itr2( itr1);
        ++itr2;
        for( ; itr2 != itr_end; ++itr2)
        {
          const t_DataType current_deviation( ( itr2->second - itr1->second).Norm());
          if( !math::EqualWithinTolerance( t_DataType( 0), current_deviation, t_DataType( 0.01), t_DataType( 0.015)))
          {
            deviation = current_deviation;
            equal &= false;
          }
        }
      }

      // print result if it was not equal within tolerance
      if( !equal)
      {
        BCL_MessageDbg( "deviation from 0.0: " + util::Format()( deviation) + "\n" + util::Format()( RESULTS));
      }

      // end
      return equal;
    }

    //! @brief check if all matrices are equal
    //! @param RESULTS map of operations string and the result matrix to be compare to other matrix
    //! @return true, if all matrices are equal within a small margin
    template< typename t_DataType>
    bool EqualWithinTolerance( const storage::Map< std::string, linal::Matrix< t_DataType> > &RESULTS) const
    {
      bool equal( true);

      // iterate over all result pairs
      for( typename storage::Map< std::string, linal::Matrix< t_DataType> >::const_iterator itr1( RESULTS.Begin()), itr_end( RESULTS.End()); itr1 != itr_end; ++itr1)
      {
        typename storage::Map< std::string, linal::Matrix< t_DataType> >::const_iterator itr2( itr1);
        ++itr2;
        for( ; itr2 != itr_end; ++itr2)
        {
          equal &= BCL_ExampleIndirectCheckWithinTolerance
                   (
                     itr2->second,
                     itr1->second,
                     t_DataType( 0.001),
                     itr1->first + " vs " + itr2->first
                   );
        }
      }

      // print result if it was not equal within tolerance
      if( !equal)
      {
        BCL_MessageDbg( util::Format()( RESULTS));
      }

      // end
      return equal;
    }

    //! @brief check if all values are equal
    //! @param RESULTS map of operations string and the result values to be compare to other values
    //! @return true, if all values are equal within a small margin
    template< typename t_DataType>
    bool EqualWithinTolerance( const storage::Map< std::string, t_DataType> &RESULTS) const
    {
      bool equal( true);

      // iterate over all result pairs
      for( typename storage::Map< std::string, t_DataType>::const_iterator itr1( RESULTS.Begin()), itr_end( RESULTS.End()); itr1 != itr_end; ++itr1)
      {
        typename storage::Map< std::string, t_DataType>::const_iterator itr2( itr1);
        ++itr2;
        for( ; itr2 != itr_end; ++itr2)
        {
          const t_DataType deviation( itr1->second - itr2->second);
          equal &= math::EqualWithinTolerance( t_DataType( 0), deviation, t_DataType( 0.01), t_DataType( 0.015));
        }
      }

      // print result if it was not equal within tolerance
      if( !equal)
      {
        BCL_MessageDbg( util::Format()( RESULTS));
      }

      // end
      return equal;
    }

  ////////////////
  // operations //
  ////////////////

    template< typename t_DataType>
    int TestAllOperations() const
    {

      // generate two Vector< double>
      t_DataType c1[ 8] = { 1, 1, 1, 1, 1, 1, 1,1 };
      t_DataType c2[ 8] = { 2, 2, 2, 2, 2, 2, 2, 2};
      t_DataType c3[ 8] = { 3, 3, 3, 3, 3, 3, 3, 3};

      linal::Vector< t_DataType> vc1( 8, c1);
      linal::Vector< t_DataType> vc2( 8, c2);
      linal::Vector< t_DataType> vc3( 8, c3);

      util::SiPtrVector< const linal::VectorInterface< t_DataType> > list_of_vectors1, list_of_vectors3;

      for( size_t i( 0); i < 10; ++i)
      {
        list_of_vectors3.PushBack( &vc1);
        list_of_vectors3.PushBack( &vc2);
        list_of_vectors3.PushBack( &vc3);
      }

      list_of_vectors1.PushBack( &vc1);
      list_of_vectors1.PushBack( &vc2);
      list_of_vectors1.PushBack( &vc3);
      list_of_vectors1.PushBack( &vc1);
      list_of_vectors1.PushBack( &vc2);
      list_of_vectors1.PushBack( &vc3);
      list_of_vectors1.PushBack( &vc1);
      list_of_vectors1.PushBack( &vc2);

      util::SiPtrVector< const linal::VectorInterface< t_DataType> > list_of_vectors2;
      list_of_vectors2.PushBack( &vc2);
      list_of_vectors2.PushBack( &vc3);
      list_of_vectors2.PushBack( &vc1);
      list_of_vectors2.PushBack( &vc2);
      list_of_vectors2.PushBack( &vc3);
      list_of_vectors2.PushBack( &vc1);
      list_of_vectors2.PushBack( &vc2);
      list_of_vectors2.PushBack( &vc3);

      // construct vector for sgemv
      const size_t vector_size_new( 325);
      const size_t matrix_height( 37);
      const size_t matrix_width( 99);

      // construct matrices from rows, cols
      linal::Matrix< t_DataType> mat1( matrix_height, vector_size_new, 2.0f);
      linal::Vector< t_DataType> vec1( vector_size_new, 3.0f);
      std::generate( vec1.Begin(), vec1.End(), Increment< t_DataType>());
      linal::Vector< t_DataType> vec2( vector_size_new, 3.0f);
      linal::Matrix< t_DataType> mat2( vector_size_new, matrix_width, 2.0f);
      std::generate( mat2.Begin(), mat2.End(), Increment< t_DataType>());

    ////////////////
    // operations //
    ////////////////

      // table that stores the timings
      storage::Table< std::string> timing_table( storage::TableHeader( linal::GetOperations< t_DataType>().GetEnumStrings()));

      {
        storage::Map< std::string, t_DataType> results;
        storage::Row< std::string> &dot_row( timing_table.InsertRow( "dot product"));
        for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
        {
          util::Stopwatch stopwatch( op_itr->GetName());
          const t_DataType dot_product( ( **op_itr)->DotProduct( vec1, vec2));
          results[ op_itr->GetName()] = dot_product;
          BCL_MessageStd( op_itr->GetName() + " dot product: " + util::Format()( dot_product));
          dot_row[ op_itr->GetName()] = stopwatch.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds();
        }
        BCL_ExampleIndirectCheck( EqualWithinTolerance< t_DataType>( results), true, "dot product results equal within tolerance");
      }

      {
        storage::Map< std::string, t_DataType> results;
        storage::Row< std::string> &norm_row( timing_table.InsertRow( "vector norm"));
        for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
        {
          util::Stopwatch stopwatch( op_itr->GetName());
          const t_DataType norm( ( **op_itr)->Norm( vec1));
          results[ op_itr->GetName()] = norm;
          BCL_MessageStd( op_itr->GetName() + " norm: " + util::Format()( norm));
          norm_row[ op_itr->GetName()] = stopwatch.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds();
        }
        BCL_ExampleIndirectCheck( EqualWithinTolerance< t_DataType>( results), true, "norm results equal within tolerance");
      }

      {
        storage::Map< std::string, linal::Vector< t_DataType> > results;
        storage::Row< std::string> &sgemv_row( timing_table.InsertRow( " sgemv [" + util::Format()( mat1.GetNumberRows()) + "," + util::Format()( mat1.GetNumberCols()) + "] x [" + util::Format()( vec1.GetSize()) + "]"));
        for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
        {
          util::Stopwatch stopwatch( op_itr->GetName());
          const linal::Vector< t_DataType> sgemv( ( **op_itr)->Multiply( mat1, vec1));
          results[ op_itr->GetName()] = sgemv;
          BCL_MessageStd( op_itr->GetName() + " sgemv: " + util::Format()( sgemv));
          sgemv_row[ op_itr->GetName()] = stopwatch.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds();
        }
        BCL_ExampleIndirectCheck( EqualWithinTolerance< t_DataType>( results), true, "sgemv results equal within tolerance");
      }

      {
        storage::Map< std::string, linal::Matrix< t_DataType> > results;
        storage::Row< std::string> &sgemm_row( timing_table.InsertRow( " sgemm [" + util::Format()( mat1.GetNumberRows()) + "," + util::Format()( mat1.GetNumberCols()) + "] x [" + util::Format()( mat2.GetNumberRows()) + "," + util::Format()( mat2.GetNumberCols()) + "]"));
        for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
        {
          util::Stopwatch stopwatch( op_itr->GetName());
          const linal::Matrix< t_DataType> sgemm( ( **op_itr)->Multiply( mat1, mat2));
          results[ op_itr->GetName()] = sgemm;
          sgemm_row[ op_itr->GetName()] = stopwatch.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds();
        }
        BCL_ExampleIndirectCheck( EqualWithinTolerance< t_DataType>( results), true, "sgemm results equal within tolerance");
      }

      {
        storage::Map< std::string, t_DataType> results;
        storage::Row< std::string> &sum_row( timing_table.InsertRow( " sum [" + util::Format()( vec1.GetSize()) + "]"));
        for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
        {
          util::Stopwatch stopwatch( op_itr->GetName());
          const t_DataType sum( ( **op_itr)->Sum( vec1));
          results[ op_itr->GetName()] = sum;
          sum_row[ op_itr->GetName()] = stopwatch.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds();
        }
        BCL_ExampleIndirectCheck( EqualWithinTolerance< t_DataType>( results), true, "sum results equal within tolerance");
      }

      {
        storage::Map< std::string, t_DataType> results;
        storage::Row< std::string> &min_row( timing_table.InsertRow( " min [" + util::Format()( vec1.GetSize()) + "]"));
        for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
        {
          util::Stopwatch stopwatch( op_itr->GetName());
          const t_DataType min( ( **op_itr)->Min( vec1));
          results[ op_itr->GetName()] = min;
          min_row[ op_itr->GetName()] = stopwatch.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds();
        }
        BCL_ExampleIndirectCheck( EqualWithinTolerance< t_DataType>( results), true, "min results equal within tolerance");
      }

      {
        storage::Map< std::string, t_DataType> results;
        storage::Row< std::string> &max_row( timing_table.InsertRow( " max [" + util::Format()( vec1.GetSize()) + "]"));
        for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
        {
          util::Stopwatch stopwatch( op_itr->GetName());
          const t_DataType max( ( **op_itr)->Max( vec1));
          results[ op_itr->GetName()] = max;
          max_row[ op_itr->GetName()] = stopwatch.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds();
        }
        BCL_ExampleIndirectCheck( EqualWithinTolerance< t_DataType>( results), true, "max results equal within tolerance");
      }

      for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
      {
        util::Stopwatch stopwatch( op_itr->GetName());
        linal::Matrix< t_DataType> square_distance_matrix( ( **op_itr)->DistanceMatrix( list_of_vectors1));
        BCL_MessageStd( op_itr->GetName() + " DistanceMatrix check for one matrix: " + util::Format()( square_distance_matrix));
      }

      for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
      {
        util::Stopwatch stopwatch( op_itr->GetName());
        linal::Matrix< t_DataType> square_distance_matrix( ( **op_itr)->DistanceMatrix( list_of_vectors3));
        BCL_MessageStd( op_itr->GetName() + " DistanceMatrix timing test for one matrix of [30][3]: ");
      }

      for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
      {
        util::Stopwatch stopwatch( op_itr->GetName());
        linal::Matrix< t_DataType> square_distance_matrix( ( **op_itr)->DistanceMatrix( list_of_vectors1, list_of_vectors2));
        BCL_MessageStd( op_itr->GetName() + " DistanceMatrix check for two matrixes: " + util::Format()( square_distance_matrix));
      }

      for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
      {
        util::Stopwatch stopwatch( op_itr->GetName());
        t_DataType result( ( **op_itr)->Sum( vec1));
        BCL_MessageStd( op_itr->GetName() + " Sum result test for one vector of 325: " + util::Format()( result));
      }

      for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
      {
        util::Stopwatch stopwatch( op_itr->GetName());
        t_DataType result( ( **op_itr)->Min( vec1));
        BCL_MessageStd( op_itr->GetName() + " Min result test for one vector of 325: " + util::Format()( result));
      }

      for( typename linal::Operations< t_DataType>::const_iterator op_itr( linal::GetOperations< t_DataType>().Begin()), op_itr_end( linal::GetOperations< t_DataType>().End()); op_itr != op_itr_end; ++op_itr)
      {
        util::Stopwatch stopwatch( op_itr->GetName());
        t_DataType result( ( **op_itr)->Max( vec1));
        BCL_MessageStd( op_itr->GetName() + " Max result test for one vector of 325: " + util::Format()( result));
      }

      // write timings
      BCL_MessageStd( "timings:");
      timing_table.WriteFormatted( util::GetLogger());

      BCL_MessageStd( "Using operations for float: " + linal::GetOperations< float>().GetDefaultOperations().GetClassIdentifier());
      BCL_MessageStd( "Using operations for double: " + linal::GetOperations< double>().GetDefaultOperations().GetClassIdentifier());

    ///////////////////////
    // general functions //
    ///////////////////////

      return 0;
    }

    int Run() const
    {
      int error( 0);
      error += TestAllOperations< float>();
      error += TestAllOperations< double>();

      return error;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalOperations

  const ExampleClass::EnumType ExampleLinalOperations::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalOperations())
  );

} // namespace bcl
