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
#include "descriptor/bcl_descriptor_limit.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "io/bcl_io_ifstream.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_limit.cpp
  //!
  //! @author mendenjl
  //! @date Oct 22, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorLimit :
    public ExampleInterface
  {
  public:

    ExampleDescriptorLimit *Clone() const
    {
      return new ExampleDescriptorLimit( *this);
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

      // test default constructor
      descriptor::Limit< chemistry::AtomConformationalInterface> small_mol_prop_log_default;

      // test Clone
      util::ShPtr< descriptor::Limit< chemistry::AtomConformationalInterface> > sp_small_mol_prop_log( small_mol_prop_log_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + sp_small_mol_prop_log->GetClassIdentifier());
      BCL_ExampleCheck
      (
        GetStaticClassName< descriptor::Limit< chemistry::AtomConformationalInterface> >(),
        small_mol_prop_log_default.GetClassIdentifier()
      );

      // Need to set up code
      BCL_MessageDbg( "begin setting up code");

      //Initialize QSAR code and set desired properties.
      descriptor::Combine< chemistry::AtomConformationalInterface, float> qsar_code;

      qsar_code.PushBack( descriptor::CheminfoProperty( util::ObjectDataLabel( "Atom_Mass")));
      qsar_code.PushBack( descriptor::CheminfoProperty( util::ObjectDataLabel( "Constant(5)")));
      qsar_code.PushBack( descriptor::CheminfoProperty( util::ObjectDataLabel( "Constant(0)")));

      BCL_MessageDbg( "Finished QSAR code setup");

      // set member variable m_LogProperty
      small_mol_prop_log_default.SetProperty( qsar_code);

      BCL_MessageDbg( "Log code: " + util::Format()( small_mol_prop_log_default));

    ////////////////
    // operations //
    ////////////////

      // setup input stream
      io::IFStream input_sdf;

      // read in molecule
      const std::string filename_in( AddExampleInputPathToFilename( e_Chemistry, "benzene_prepared.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, filename_in);
      chemistry::FragmentComplete small_mol = sdf::FragmentFactory::MakeFragment( input_sdf);
      BCL_MessageStd
      (
        "pre code: " + small_mol_prop_log_default.GetString()
      );

      // test the actual code generation
      BCL_MessageStd
      (
        "The code: " + util::Format()( small_mol_prop_log_default.CollectValuesOnEachElementOfObject( small_mol))
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test the stream operators
      BCL_MessageStd
      (
        "Outputting sp_small_mol_prop_log: " + util::Format()( small_mol_prop_log_default)
      );

      // Write object file for logcode
      WriteBCLObject( small_mol_prop_log_default);

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorLimit

  const ExampleClass::EnumType ExampleDescriptorLimit::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorLimit())
  );

} // namespace bcl
