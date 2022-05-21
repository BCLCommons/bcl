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

#ifndef BCL_SCORE_PROTEIN_MODEL_AA_NEIGHBORHOOD_H_
#define BCL_SCORE_PROTEIN_MODEL_AA_NEIGHBORHOOD_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_neighborhood_interface.h"
#include "bcl_score_protein_model.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelAANeighborhood
    //! @brief a class that scores the neighborhood of each amino acid within a protein model
    //! @details This class collects all amino acids in the protein model and calls a
    //! AANeighborhoodInterface derived class with the amino acid neighbor list, the amino acids in the chain and all
    //! amino acids in different chains in order to evaluate the neighborhood potential of that very amino acid. These
    //! is summed up over all amino acids in the given protein model.
    //!
    //! @see @link example_score_protein_model_aa_neighborhood.cpp @endlink
    //! @author woetzen
    //! @date 16.01.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelAANeighborhood :
      public ProteinModel
    {

    public:

      //! normalization type
      enum NormalizationType
      {
        e_None,      // don't normalize
        e_Normalize, // divide sum by number of summed residues
        e_RMSD,      // calculate an RMSD instead of an average
        e_FractionExplained, // Amount of the prediction that is explained / contained in the model
        e_CorrelationPlusRelativeError, // Spearman correlation and relative error
        s_NumberTypes
      };

      //! @brief conversion to a string from a Type
      //! @param TYPE the type to get a string for
      //! @return a string representing that type
      static const std::string &GetTypeName( const NormalizationType &TYPE);

      //! @brief conversion to a string from a Type
      //! @param TYPE the type to get a string for
      //! @return a string representing that type
      static const std::string &GetTypeSuffix( const NormalizationType &TYPE);

      //! @brief enum class wrapper for Type
      typedef util::WrapperEnum< NormalizationType, &GetTypeName, s_NumberTypes> TypeEnum;

    private:

    //////////
    // data //
    //////////

      //! AAExposure function to be used for calculations
      util::ShPtr< AANeighborhoodInterface> m_ScoreAANeighborhood;

      //! ShPtr to AANeighborListContainer generator
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > m_AANeighborListContainerGenerator;

      //! normalization type
      TypeEnum m_Normalization;

      //! scheme to be used
      std::string m_Scheme;

      //! score type
      ProteinModel::Type m_ScoreType;

      //! readable scheme
      std::string m_ReadableScheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelAANeighborhood();

      //! @brief constructor from ShPtr to aa exposure function
      //! @param SP_SCORE_AA_NEIGHBORHOOD ShPtr to AA neighborhood scoring function
      //! @param TYPE normalization type to use
      //! @param CONSIDER_DIFFERENT_CHAIN whether neighbors from different chains should be considered
      //! @param SCORE_TYPE score type
      //! @param READABLE_SCHEME scheme that is more human readable
      ProteinModelAANeighborhood
      (
        const util::ShPtr< AANeighborhoodInterface> &SP_SCORE_AA_NEIGHBORHOOD,
        const TypeEnum &TYPE = e_None,
        const bool CONSIDER_DIFFERENT_CHAIN = true,
        const ProteinModel::Type &SCORE_TYPE = ProteinModel::e_Undefined,
        const std::string &READABLE_SCHEME = ""
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new ProteinModelAANeighborhood copied from this one
      ProteinModelAANeighborhood *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      const std::string &GetReadableScheme() const
      {
        return m_ReadableScheme;
      }

      //! @brief get score type
      //! @return score type
      Type GetType() const
      {
        return m_ScoreType;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the sum of exposures of all amino acids for the given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return the sum of exposures of all amino acids for the given ProteinModel
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM
      ) const;

    }; // class ProteinModelAANeighborhood

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_PROTEIN_MODEL_AA_NEIGHBORHOOD_H_
