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
#include "model/bcl_model_descriptor_selection_backward_elimination.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_descriptor_selection_backward_elimination.cpp
  //!
  //! @author butkiem1
  //! @date May 31, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDescriptorSelectionBackwardElimination :
    public ExampleInterface
  {
  public:

    ExampleModelDescriptorSelectionBackwardElimination *Clone() const
    {
      return new ExampleModelDescriptorSelectionBackwardElimination( *this);
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
      // sample empty descriptor set
      util::ObjectDataLabel empty_descriptor_set( "Combine");

      // sample initial descriptor set
      util::ObjectDataLabel initial_descriptor_set
      (
        "Combine("
        "  Sum(Atom_Mass),"
        "  Sum(Atom_HbondDonors),"
        "  Sum(Atom_HbondAcceptors)"
        ")"
      );

      // sample entire descriptor set
      util::ObjectDataLabel entire_descriptor_set
      (
        "Combine("
        "  Sum(Atom_Mass),"
        "  Sum(Atom_HbondDonors),"
        "  Sum(Atom_HbondAcceptors),"
        "  LogP"
        ")"
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::DescriptorSelectionBackwardElimination selection_default;

      // construct through object label
      // initialize descriptor selection scheme
      util::Implementation< model::DescriptorSelectionInterface> descriptor_selection
      (
        util::ObjectDataLabel( std::string( "FeatureBackwardElimination"))
      );

    /////////////////
    // data access //
    /////////////////

      // check alias
      BCL_ExampleCheck( selection_default.GetAlias(), std::string( "FeatureBackwardElimination"));
      BCL_ExampleCheck( descriptor_selection.GetAlias(), std::string( "FeatureBackwardElimination"));

      // check GetInitialDescriptorSet
      BCL_ExampleCheck( selection_default.GetInitialDescriptorSet( entire_descriptor_set), entire_descriptor_set);
      BCL_ExampleCheck( descriptor_selection->GetInitialDescriptorSet( entire_descriptor_set), entire_descriptor_set);

    ///////////////
    // operators //
    ///////////////

      BCL_ExampleCheck
      (
        selection_default( empty_descriptor_set, entire_descriptor_set).GetSize(), size_t( 0)
      );

      BCL_ExampleCheck
      (
        selection_default( initial_descriptor_set, entire_descriptor_set).GetSize(), size_t( 3)
      );

      BCL_ExampleCheck
      (
        descriptor_selection->operator ()( empty_descriptor_set, entire_descriptor_set).GetSize(), size_t( 0)
      );

      BCL_ExampleCheck
      (
        descriptor_selection->operator ()( initial_descriptor_set, entire_descriptor_set).GetSize(), size_t( 3)
      );

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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDescriptorSelectionBackwardElimination

  const ExampleClass::EnumType ExampleModelDescriptorSelectionBackwardElimination::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDescriptorSelectionBackwardElimination())
  );

} // namespace bcl
