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
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_constitution_graph_converter.cpp
  //!
  //! @author kothiwsk
  //! @date Jan 24, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConstitutionGraphConverter :
    public ExampleInterface
  {
  public:

    ExampleChemistryConstitutionGraphConverter *Clone() const
    {
      return new ExampleChemistryConstitutionGraphConverter( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      chemistry::ConstitutionGraphConverter default_constructor;

      // element type constructor
      chemistry::ConstitutionGraphConverter constructor_a
      (
        chemistry::ConstitutionGraphConverter::e_ElementType,
        chemistry::ConstitutionalBondTypeData::e_BondOrderOrAromaticWithRingness
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // read in some molecules
      io::IFStream input_sdf;
      const std::string taxol_filename( AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, taxol_filename);

      // load information into small_mol_conformation
      chemistry::FragmentComplete small_mol_conformation;
      small_mol_conformation
        = sdf::FragmentFactory::MakeFragment( input_sdf);

      // close the input stream
      io::File::CloseClearFStream( input_sdf);

      // create constitution
      chemistry::FragmentConstitutionShared constitution( small_mol_conformation);

      // create graphs using constitution
      graph::ConstGraph< size_t, size_t> const_graph_taxol( default_constructor( constitution));
      graph::ConstGraph< size_t, size_t> const_graph_taxol_a( constructor_a( constitution));

      BCL_ExampleCheck( const_graph_taxol.GetSize(), constitution.GetNumberAtoms());
      BCL_ExampleCheck( const_graph_taxol.NumEdges(), 2 * constitution.GetNumberBonds());
      BCL_ExampleCheck( const_graph_taxol.AreConnected( 0, 1), true);
      BCL_ExampleCheck( !const_graph_taxol.AreConnected( 0, 6), true);
      BCL_ExampleCheck( graph::Connectivity::IsConnected( const_graph_taxol), true);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryConstitutionGraphConverter

  const ExampleClass::EnumType ExampleChemistryConstitutionGraphConverter::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConstitutionGraphConverter())
  );

} // namespace bcl
