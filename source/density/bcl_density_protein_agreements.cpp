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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "density/bcl_density_protein_agreements.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom_types.h"
#include "density/bcl_density_protein_agreement_ccc.h"
#include "density/bcl_density_protein_agreement_likelihood.h"
#include "density/bcl_density_simulators.h"
#include "math/bcl_math_polynomial.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    const double g_LikelihoodCaMeanPolynomCoeff[]       = {  0.00000346634, -0.000149827,   0.00159432, 0.00788547, 0.00386209};
    const double g_LikelihoodCaSdPolynomCoeff[]         = { -0.000000465365, 0.0000547865, -0.00153436, 0.0161681,  0.0339695};
    const double g_LikelihoodCbMeanPolynomCoeff[]       = {  0.00000322947, -0.000139397,   0.00145854, 0.00815733, 0.00711573};
    const double g_LikelihoodCbSdPolynomCoeff[]         = { -0.00000126792,  0.0000982248, -0.00234045, 0.0217593,  0.0257755};
    const double g_LikelihoodBackBoneMeanPolynomCoeff[] = {  0.00000427199,  0.000235853,  -0.00546931, 0.0659212, -0.12085};
    const double g_LikelihoodBackBoneSdPolynomCoeff[]   = { -0.00005495,     0.000288613,  -0.00544743, 0.044201,  -0.021971};

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinAgreements::ProteinAgreements() :
      e_CCC(                AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementCCC( false, false)))),
      e_CCCScaled(          AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementCCC( false, true)))),
      e_LikelihoodCa(       AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementLikelihood( false, storage::Set< biol::AtomType>( biol::GetAtomTypes().CA), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCaMeanPolynomCoeff)), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCaSdPolynomCoeff)))))),
      e_LikelihoodCb(       AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementLikelihood( false, storage::Set< biol::AtomType>( biol::GetAtomTypes().CB), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCbMeanPolynomCoeff)), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodCbSdPolynomCoeff)))))),
      e_LikelihoodBackBone( AddAgreement( util::ShPtr< ProteinAgreementInterface>( new ProteinAgreementLikelihood( false, biol::GetAtomTypes().GetBackBoneAtomTypes()            , math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodBackBoneMeanPolynomCoeff)), math::Polynomial::MakeFromCoefficients( linal::Vector< double>( 5, g_LikelihoodBackBoneSdPolynomCoeff))))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinAgreements::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a new protein agreement instance
    //! @param SP_PROTEIN_AGREEMENT shptr to instance of a ProteinAgreementInterface derived class
    //! @return the new enum constructed from the agreement
    ProteinAgreement ProteinAgreements::AddAgreement( const util::ShPtr< ProteinAgreementInterface> &SP_PROTEIN_AGREEMENT)
    {
      return AddEnum( SP_PROTEIN_AGREEMENT->GetScheme(), SP_PROTEIN_AGREEMENT);
    }

    //! @brief create ProteinAgreement from enum, density map and resolution
    //! @param PROTEIN_AGREEMENT enum of the agreement to be used
    //! @param SIMULATOR
    //! @param SP_DENSITY density map to be used
    //! @param RESOLUTION resolution in Angstrom to be used
    //! @return ShPtr to density::ProteinAgreement
    util::ShPtr< ProteinAgreementInterface> ProteinAgreements::CreateProteinAgreement
    (
      const ProteinAgreement &PROTEIN_AGREEMENT,
      const Simulator &SIMULATOR,
      const util::SiPtr< const Map> &SP_DENSITY,
      const double RESOLUTION
    ) const
    {
      // check that the enum is defined
      if( !PROTEIN_AGREEMENT.IsDefined() || !PROTEIN_AGREEMENT->IsDefined())
      {
        return util::ShPtr< ProteinAgreementInterface>();
      }

      // create simulator
      util::ShPtr< SimulateInterface> sp_simulator( GetSimulators().CreateSimulator( SIMULATOR, SP_DENSITY->GetCellWidth(), RESOLUTION));

      // clone an agreement
      util::ShPtr< ProteinAgreementInterface> sp_protein_agreement( PROTEIN_AGREEMENT->HardCopy());

      // set the member
      sp_protein_agreement->SetDensityMap( SP_DENSITY);
      sp_protein_agreement->SetSimulator( sp_simulator);

      // end
      return sp_protein_agreement;
    }

    //! @brief construct on access function for all ProteinAgreements
    //! @return reference to only instances of ProteinAgreements
    ProteinAgreements &GetProteinAgreements()
    {
      return ProteinAgreements::GetEnums();
    }

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< density::ProteinAgreementInterface>, density::ProteinAgreements>;

  } // namespace util
} // namespace bcl
