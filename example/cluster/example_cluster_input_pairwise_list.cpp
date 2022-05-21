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
#include "cluster/bcl_cluster_input_pairwise_list.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    // explicit instantiation of List< std::string>
    template class List< std::string>;
  } // namespace storage

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_input_pairwise_list.cpp
  //!
  //! @author alexanns
  //! @date Nov. 30, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterInputPairwiseList :
    public ExampleInterface
  {
  public:

      ExampleClusterInputPairwiseList *Clone() const
    {
      return new ExampleClusterInputPairwiseList( *this);
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

    int Run() const
    {
      // set up some data
      util::ShPtr< storage::List< std::string> > object_list
      (
        new storage::List< std::string>()
      );
      std::string one( "one");
      std::string two( "two");
      std::string three( "three");

      object_list->PushBack( one);
      object_list->PushBack( two);
      object_list->PushBack( three);

      storage::HashMap< size_t, storage::HashMap< size_t, double> > data;
      data[ 2][ 1] = 21.0;
      data[ 1][ 3] = 13.0;
      data[ 2][ 3] = 23;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // check default constructor
      cluster::InputPairwiseList< double> default_constr;

      const std::string object_list_filename
      (
        AddExampleInputPathToFilename( e_Cluster, "input_pairwise_list_object_list.txt")
      );

      // check that the object was created properly by handling an input file
      io::IFStream read;
      const std::string input_pairwise_filename
      (
        AddExampleInputPathToFilename( e_Cluster, "input_pairwise_list.ls")
      );
      BCL_ExampleMustOpenInputFile( read, input_pairwise_filename);
      util::ShPtr< cluster::DistancesStored< std::string, double> > handled_output( default_constr.HandleInput( read));
      io::File::CloseClearFStream( read);

      // make sure the size of the containers are correct
      BCL_ExampleCheck( default_constr.GetInputObjects()->GetSize(), object_list->GetSize());
      BCL_ExampleCheck( handled_output->GetData().GetSize(), data.GetSize());

      // make sure that the read in objects are correct
      {
        storage::List< std::string>::const_iterator object_a_itr
        (
          std::find( default_constr.GetInputObjects()->Begin(), default_constr.GetInputObjects()->End(), two)
        );
        storage::List< std::string>::const_iterator object_b_itr
        (
          std::find( default_constr.GetInputObjects()->Begin(), default_constr.GetInputObjects()->End(), one)
        );
        storage::HashMap< size_t, storage::HashMap< size_t, double> > input_data( handled_output->GetData());
        const double pairwise_value( input_data[ size_t( &( *object_a_itr))][ size_t( &( *object_b_itr))]);
        const double correct_value( 21);
        BCL_ExampleCheck( pairwise_value, correct_value);
      }
      {
        storage::List< std::string>::const_iterator object_a_itr
        (
          std::find( default_constr.GetInputObjects()->Begin(), default_constr.GetInputObjects()->End(), one)
        );
        storage::List< std::string>::const_iterator object_b_itr
        (
          std::find( default_constr.GetInputObjects()->Begin(), default_constr.GetInputObjects()->End(), three)
        );
        storage::HashMap< size_t, storage::HashMap< size_t, double> > input_data( handled_output->GetData());
        const double pairwise_value( input_data[ size_t( &( *object_a_itr))][ size_t( &( *object_b_itr))]);
        const double correct_value( 13);
        BCL_ExampleCheck( pairwise_value, correct_value);
      }
      {
        storage::List< std::string>::const_iterator object_a_itr
        (
          std::find( default_constr.GetInputObjects()->Begin(), default_constr.GetInputObjects()->End(), two)
        );
        storage::List< std::string>::const_iterator object_b_itr
        (
          std::find( default_constr.GetInputObjects()->Begin(), default_constr.GetInputObjects()->End(), three)
        );
        storage::HashMap< size_t, storage::HashMap< size_t, double> > input_data( handled_output->GetData());
        const double pairwise_value( input_data[ size_t( &( *object_a_itr))][ size_t( &( *object_b_itr))]);
        const double correct_value( 23);
        BCL_ExampleCheck( pairwise_value, correct_value);
      }

      // test copy constructor
      cluster::InputPairwiseList< double> copy_constr( default_constr);

      // test copy constructor by handling the input from the file
      BCL_ExampleMustOpenInputFile( read, input_pairwise_filename);

      util::ShPtr< cluster::DistancesStored< std::string, double> > handled_output_copy_constr( copy_constr.HandleInput( read));
      io::File::CloseClearFStream( read);

      // make sure the size of the containers are correct
      BCL_ExampleCheck( copy_constr.GetInputObjects()->GetSize(), object_list->GetSize());
      BCL_ExampleCheck( handled_output_copy_constr->GetData().GetSize(), data.GetSize());

      // test clone constructor
      util::ShPtr< cluster::InputPairwiseList< double> > clone_constr( default_constr.Clone());

      // test copy constructor by handling the input from the file
      BCL_ExampleMustOpenInputFile( read, input_pairwise_filename);
      util::ShPtr< cluster::DistancesStored< std::string, double> > handled_output_clone_constr
      (
        clone_constr->HandleInput( read)
      );
      io::File::CloseClearFStream( read);

      // make sure the containers from the clone constructor have the correct size
      BCL_ExampleCheck( clone_constr->GetInputObjects()->GetSize(), object_list->GetSize());
      BCL_ExampleCheck( handled_output_clone_constr->GetData().GetSize(), data.GetSize());

    /////////////////
    // data access //
    /////////////////

      const std::string correct_static_class_name( "bcl::cluster::InputPairwiseList<double>");
      // check GetClassIdentifier
      BCL_ExampleCheck( correct_static_class_name, clone_constr->GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // HandleInput function is tested and demonstrated while checking constructors

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      WriteBCLObject( default_constr);

      // read the file back
      BCL_MessageStd( "Testing the read function");
      cluster::InputPairwiseList< double> read_output;
      ReadBCLObject( read_output);
      BCL_MessageDbg( "read in object is  \n" + util::Format()( read_output));

      // make sure the read function worked correctly by handling an input file with the read in object
      BCL_ExampleMustOpenInputFile( read, input_pairwise_filename);
      util::ShPtr< cluster::DistancesStored< std::string, double> > handled_output_read( read_output.HandleInput( read));
      io::File::CloseClearFStream( read);
      BCL_MessageDbg( "handled output from written object is \n" + util::Format()( handled_output));
      BCL_MessageDbg
      (
        "handled output from read in object is \n" + util::Format()( handled_output_read)
      );

      // make sure the read in object has the correct size of containers
      BCL_ExampleCheck( read_output.GetInputObjects()->GetSize(), object_list->GetSize());
      BCL_ExampleCheck( handled_output_read->GetData().GetSize(), data.GetSize());

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterInputPairwiseList

  const ExampleClass::EnumType ExampleClusterInputPairwiseList::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterInputPairwiseList())
  );
} // namespace bcl
