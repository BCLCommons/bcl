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
#include "model/bcl_model_kohonen_node.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_kohonen_node.cpp
  //!
  //! @author lemmonwa, mueller
  //! @date Aug 26, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelKohonenNode :
    public ExampleInterface
  {
  public:

    ExampleModelKohonenNode *Clone() const
    {
      return new ExampleModelKohonenNode( *this);
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
      const float feature_array[] = { 34.5, 20.1, 42.0};
      const float result_array[] = { 1.0};
      const float position_array[] = { 1.0, 2.0};

      // set up parts for  testing
      linal::Vector< float> feature( 3, feature_array), result( 1, result_array), position( 2, position_array);
      linal::Vector< float> zero_feature( 3, 0.0), zero_result( 1, 0.0), zero_position( 2, 0.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::KohonenNode default_node;

      // constructor from all parameters
      model::KohonenNode node( position, feature, result);

      // construct a node initialized with a zero reference vector
      model::KohonenNode zero_node( zero_position, zero_feature, zero_result);

    /////////////////
    // data access //
    /////////////////

      // test inline linal::Vector< float> &GetFeatureVector()
      BCL_ExampleCheck( node.GetFeatureVector(), feature);

      // test inline linal::Vector< float> &GetResultVector()
      BCL_ExampleCheck( node.GetResultVector(), result);

      // test inline const linal::Vector< size_t> &GetPosition() const
      BCL_ExampleCheck( node.GetPosition(), position);

      // test GetMappedData() const
      BCL_ExampleCheck( node.GetWeight(), 0.0);

    ////////////////
    // operations //
    ////////////////

      // test void InitializeMap()
      node.MapData( feature, result);
      node.Reset();
      BCL_ExampleCheck( node.GetWeight(), 0.0);

      // test void MapData( const storage::VectorND< 2, linal::Vector< float> > &DATA)
      node.MapData( feature, result);
      BCL_ExampleCheck( node.GetFeatureVector(), feature);
      BCL_ExampleCheck( node.GetResultVector(), result);
      BCL_ExampleCheck( node.GetWeight(), double( 1.0));

    ///////////////
    // operators //
    ///////////////

      // test operator+=( const KohonenNode &OTHER);
      node += node;
      BCL_ExampleCheck( node.GetFeatureVector(), feature);
      BCL_ExampleCheck( node.GetResultVector(), result);
      BCL_ExampleCheck( node.GetWeight(), double( 2.0));

      // test operator=( const KohonenNode &OTHER);
      node = zero_node;
      BCL_ExampleCheck( node.GetWeight(), 0.0);

    //////////////////////
    // input and output //
    //////////////////////

      // write bcl object
      WriteBCLObject( node);

      // create default object
      model::KohonenNode node_read;

      // read bcl object
      ReadBCLObject( node_read);

      // check read in object
      BCL_ExampleCheck( node_read.GetPosition(), node.GetPosition());

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleClass::ExampleModelKohonenNode

  const ExampleClass::EnumType ExampleModelKohonenNode::s_Instance
  (
    GetExamples().AddEnum( ExampleModelKohonenNode())
  );

} // namespace bcl
