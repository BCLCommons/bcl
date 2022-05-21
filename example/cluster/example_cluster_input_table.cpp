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
#include "cluster/bcl_cluster_input_table.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_input_table.cpp
  //!
  //! @author alexanns
  //! @date Sep 09, 2009
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterInputTable :
    public ExampleInterface
  {
  public:

    ExampleClusterInputTable *Clone() const
    {
      return new ExampleClusterInputTable( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "Testing default constructor");
      cluster::InputTable< double> default_const;

      // constructor indicating whole table must be read in
      BCL_MessageStd( "constructor indicating whole table must be read in");
      cluster::InputTable< double> entire_table_const( true, true);

      // constructor indicating upper triangle of table must be read in
      BCL_MessageStd( "constructor indicating upper triangle of table must be read in");
      cluster::InputTable< double> upper_triangle_table_const( true, false);

      // constructor indicating lower triangle of table must be read in
      BCL_MessageStd( "constructor indicating lower triangle of table must be read in");
      cluster::InputTable< double> lower_triangle_table_const( false, true);

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // test reading entire matrix :
      {
        // create string "table_entire_filename" and initialize with the filename of the entire matrix
        std::string table_filename( AddExampleInputPathToFilename( e_Cluster, "table_entire.txt"));

        // create ifstream read
        io::IFStream read;

        // open "table_entire_filename"
        BCL_ExampleMustOpenInputFile( read, table_filename);

        // create "member_distance_function" which will be used to calculate distances between members of a node
        // so that the center of the node can be calculated
        util::ShPtr< cluster::DistancesStored< std::string, double> > handled_input( entire_table_const.HandleInput( read));

        io::File::CloseClearFStream( read);

        {
          const size_t expected_size( 4);
          const size_t actual_size( handled_input->GetData().GetSize());
          BCL_ExampleCheck( actual_size, expected_size);
        }
        {
          const size_t expected_size( 8);
          const size_t actual_size( entire_table_const.GetInputObjects()->GetSize());
          BCL_ExampleCheck( actual_size, expected_size);
        }
        {
          const size_t expected_value( 1);
          const size_t search_a( size_t( &( *++++++++entire_table_const.GetInputObjects()->Begin())));
          const size_t search_b( size_t( &( *entire_table_const.GetInputObjects()->Begin())));
          storage::HashMap< size_t, storage::HashMap< size_t, double> >::const_iterator find_itr_a
          (
            handled_input->GetData().Find( search_a)
          );
          BCL_ExampleAssert( find_itr_a != handled_input->GetData().End(), true);
          storage::HashMap< size_t, double>::const_iterator find_itr_b
          (
            find_itr_a->second.Find( search_b)
          );
          BCL_ExampleAssert( find_itr_b != find_itr_a->second.End(), true);
          const size_t actual_value( find_itr_b->second);
          BCL_ExampleCheck( expected_value, actual_value);
        }

      }

      // test reading upper triangle matrix :
      {
        // create string "table_entire_filename" and initialize with the filename of the entire matrix
        std::string table_filename( AddExampleInputPathToFilename( e_Cluster, "table_upper_triangle.txt"));

        // create ifstream read
        io::IFStream read;

        // open "table_entire_filename"
        BCL_ExampleMustOpenInputFile( read, table_filename);

        util::ShPtr< cluster::DistancesStored< std::string, double> > handled_input( upper_triangle_table_const.HandleInput( read));

        io::File::CloseClearFStream( read);

        {
          const size_t expected_size( 3);
          const size_t actual_size( handled_input->GetData().GetSize());
          BCL_ExampleCheck( actual_size, expected_size);
        }
        {
          const size_t expected_size( 4);
          const size_t actual_size( upper_triangle_table_const.GetInputObjects()->GetSize());
          BCL_ExampleCheck( actual_size, expected_size);
        }
        {
          const size_t expected_value( 7);

          const size_t search_a( size_t( &( *++upper_triangle_table_const.GetInputObjects()->Begin())));
          const size_t search_b( size_t( &( *++++upper_triangle_table_const.GetInputObjects()->Begin())));
          storage::HashMap< size_t, storage::HashMap< size_t, double> >::const_iterator find_itr_a
          (
            handled_input->GetData().Find( search_a)
          );
          BCL_ExampleAssert( find_itr_a != handled_input->GetData().End(), true);
          storage::HashMap< size_t, double>::const_iterator find_itr_b
          (
            find_itr_a->second.Find( search_b)
          );
          BCL_ExampleAssert( find_itr_b != find_itr_a->second.End(), true);
          const size_t actual_value( find_itr_b->second);
          BCL_ExampleCheck( expected_value, actual_value);
        }
      }

      // test reading lower triangle matrix :
      {
        // create string "table_entire_filename" and initialize with the filename of the entire matrix
        std::string table_filename( AddExampleInputPathToFilename( e_Cluster, "table_lower_triangle.txt"));

        // create ifstream read
        io::IFStream read;

        // open "table_entire_filename"
        BCL_ExampleMustOpenInputFile( read, table_filename);

        util::ShPtr< cluster::DistancesStored< std::string, double> > handled_input( lower_triangle_table_const.HandleInput( read));

        io::File::CloseClearFStream( read);

        {
          const size_t expected_size( 3);
          const size_t actual_size( handled_input->GetData().GetSize());
          BCL_ExampleCheck( actual_size, expected_size);
        }
        {
          const size_t expected_size( 4);
          const size_t actual_size( lower_triangle_table_const.GetInputObjects()->GetSize());
          BCL_ExampleCheck( actual_size, expected_size);
        }
        {
          const size_t expected_value( 13);
          const size_t search_a( size_t( &( *++++++lower_triangle_table_const.GetInputObjects()->Begin())));
          const size_t search_b( size_t( &( *lower_triangle_table_const.GetInputObjects()->Begin())));
          storage::HashMap< size_t, storage::HashMap< size_t, double> >::const_iterator find_itr_a
          (
            handled_input->GetData().Find( search_a)
          );
          BCL_ExampleAssert( find_itr_a != handled_input->GetData().End(), true);
          storage::HashMap< size_t, double>::const_iterator find_itr_b
          (
            find_itr_a->second.Find( search_b)
          );
          BCL_ExampleAssert( find_itr_b != find_itr_a->second.End(), true);
          const size_t actual_value( find_itr_b->second);
          BCL_ExampleCheck( expected_value, actual_value);
        }
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterInputTable

  const ExampleClass::EnumType ExampleClusterInputTable::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterInputTable())
  );

} // namespace bcl
