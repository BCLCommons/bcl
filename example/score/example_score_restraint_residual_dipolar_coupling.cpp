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
#include "score/bcl_score_restraint_residual_dipolar_coupling.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "nmr/bcl_nmr_residual_dipolar_coupling_least_square_deviation.h"
#include "nmr/bcl_nmr_star_rdc_handler.h"
#include "score/bcl_score_residual_dipolar_coupling_q_value.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_residual_dipolar_coupling.cpp
  //!
  //! @author alexanns
  //! @date Apr 16, 2009
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintResidualDipolarCoupling :
    public ExampleInterface
  {
  public:

    ExampleScoreRestraintResidualDipolarCoupling *Clone() const
    {
      return new ExampleScoreRestraintResidualDipolarCoupling( *this);
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
      // read in the restraints
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1ubi_rdc_baxbicHN.exp"));

      // read in the protein model
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

      // get the aa sequence
      const util::SiPtrVector< const biol::AASequence> aa_seq( protein_model.GetSequences());

      // initialize the handler
      nmr::StarRDCHandler handler;
      restraint::RDC rdcs( handler.ReadRestraints( read));

      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::RestraintResidualDipolarCoupling default_constr;

      // constructor taking parameters
      score::RestraintResidualDipolarCoupling rdc_score
      (
        util::CloneToShPtr( rdcs),
        ( nmr::ResidualDipolarCouplingLeastSquareDeviation()),
        ( score::ResidualDipolarCouplingQValue())
      );

    /////////////////
    // data access //
    /////////////////

      // test GetScheme
      BCL_ExampleCheck( rdc_score.GetScheme(), "rdc_restraint");

    ///////////////
    // operators //
    ///////////////

      // test () operator
      BCL_ExampleCheckWithinTolerance( rdc_score( protein_model), -113.54, 0.0001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreRestraintResidualDipolarCoupling

  const ExampleClass::EnumType ExampleScoreRestraintResidualDipolarCoupling::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintResidualDipolarCoupling())
  );

} // namespace bcl
