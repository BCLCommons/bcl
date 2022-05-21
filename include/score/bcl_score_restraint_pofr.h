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

#ifndef BCL_SCORE_RESTRAINT_POFR_H_
#define BCL_SCORE_RESTRAINT_POFR_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_pofr.h"
#include "bcl_score_protein_model.h"
#include "restraint/bcl_restraint_sas_pofr_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintPofr
    //! @brief This Interface function will return an agreement score for a given protein model
    //! @details This Interface function uses restraint::SasPofr and score::PofR to compare the
    //! @details experimental and calculated SAS distance histograms for a given protein and report the deviation.
    //!
    //! @see @link example_score_restraint_pofr.cpp @endlink
    //! @author putnamdk
    //! @date July 16, 2015
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintPofr :
      public ProteinModel
    {

    private:

      //! storage for Saxs Intensity data from Protein Model
      util::Implementation< restraint::SasPofRInterface> m_Calc;

      PofR m_Score;

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
      RestraintPofr();

      //! @brief construct from member variables
      //! @param DATA_FROM_STRUCTURE_CALCULATOR function interface for calculating theoretical Saxs Curves
      //! @param DATA_AGREEMENT_CALCULATOR function interface for calculating a Saxs Curve agreement value
      //! @param SCHEME Scheme to be used
      RestraintPofr
      (
        const util::Implementation< restraint::SasPofRInterface> &POFR_IMPLEMENTATION,
        const PofR &SCORE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new RestraintPofr
      RestraintPofr *Clone() const;

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

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class RestraintPofr

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_POFR_H_
