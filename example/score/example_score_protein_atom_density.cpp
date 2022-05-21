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
#include "score/bcl_score_protein_atom_density.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_atom_density.cpp
  //!
  //! @author nobody
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreProteinAtomDensity :
    public ExampleInterface
  {
  public:

    ExampleScoreProteinAtomDensity *Clone() const
    { return new ExampleScoreProteinAtomDensity( *this);}

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
      //instantiate pdb
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      //instantiate proteinmodels of chains
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      score::ProteinAtomDensity atom_density
      (
        storage::Set< biol::SSType>( biol::GetSSTypes().HELIX.GetIterator(), biol::GetSSTypes().COIL.GetIterator()),
        storage::Set< biol::AtomType>( biol::GetAtomTypes().CA),
        linal::Vector3D( 5.0, 5.0, 5.0)
      );

      const double helix_ca_density( atom_density.CalculateAverageAtomDensity( model));

      BCL_MessageStd( "ca density in helices of 1CDA.pdb: " + util::Format()( helix_ca_density));

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreProteinAtomDensity

  const ExampleClass::EnumType ExampleScoreProteinAtomDensity::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinAtomDensity())
  );

} // namespace bcl
