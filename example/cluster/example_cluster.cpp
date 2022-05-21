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
#include "cluster/bcl_cluster.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCluster :
    public ExampleInterface
  {
  public:

    ExampleCluster *Clone() const
    {
      return new ExampleCluster( *this);
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

      // check if get namespace identifier is working
      const std::string namespace_identifier( cluster::GetNamespaceIdentifier());
      const std::string expected_namespace_identifier( "bcl::cluster");

      BCL_MessageStd
      (
        "this is the " + expected_namespace_identifier + " namespace identifier : " + namespace_identifier
      );

      BCL_Example_Check
      (
        namespace_identifier == expected_namespace_identifier,
        "namespace identifier for " + expected_namespace_identifier + " is incorrect. It should be "
        + expected_namespace_identifier
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorage

  const ExampleClass::EnumType ExampleCluster::s_Instance
  (
    GetExamples().AddEnum( ExampleCluster())
  );

} // namespace bcl
