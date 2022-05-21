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

#ifndef BCL_SCORE_RESTRAINT_SAXS_H_
#define BCL_SCORE_RESTRAINT_SAXS_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "bcl_score_sas_type.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintSaxs
    //! @brief This Interface function will return an RMSD score for a given protein model
    //! @details This Interface function uses restraint::SasDebye and score::SaxsRMSD to compare the
    //! @details experimental and calculated SAXS curves for a given protein and report the RMSD between them.
    //!
    //! @see @link example_score_restraint_saxs.cpp @endlink
    //! @author putnamdk
    //! @date May 12, 2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintSaxs :
      public ProteinModel
    {

    private:

      //! storage for Saxs Intensity data from Protein Model
      util::Implementation< restraint::SasDebyeInterface> m_Calc;

      SasType m_Score;

      //! scheme
      std::string m_Scheme;

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintSaxs();

      //! @brief construct from member variables
      //! @param DATA_FROM_STRUCTURE_CALCULATOR function interface for calculating theoretical Saxs Curves
      //! @param DATA_AGREEMENT_CALCULATOR function interface for calculating a Saxs Curve agreement value
      //! @param SCHEME Scheme to be used
      RestraintSaxs
      (
        const util::Implementation< restraint::SasDebyeInterface> &DEBYE_IMPLEMENTATION,
        const SasType &SCORE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new RestraintSaxs
      RestraintSaxs *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() which takes an ProteinModel for calculating its agreement with the Saxs Data
      //! @param PROTEIN_MODEL the ProteinModel which will be scored for agreement with the Saxs Data
      //! @return return a double which is the score of the agreement of the ProteinModel with the Saxs data
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    }; // class RestraintSaxs

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_SAXS_H_
