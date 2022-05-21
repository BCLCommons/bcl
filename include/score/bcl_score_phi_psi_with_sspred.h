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

#ifndef BCL_SCORE_PHI_PSI_WITH_SSPRED_H_
#define BCL_SCORE_PHI_PSI_WITH_SSPRED_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "sspred/bcl_sspred.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_phi_psi.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PhiPsiWithSSPred
    //! @brief scores the phi and psi angles of an sse taking into account secondary structure prediction methods
    //! @details The score is independent of the type of SSE. It determines the score of an amino acid according to
    //!          the best out of the phi psi score for every sstype weighted by the ss prediction for that type.
    //!
    //! @see @link example_score_phi_psi_with_sspred.cpp @endlink
    //! @author alexanns
    //! @date Jan 6, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PhiPsiWithSSPred :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> >
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! map of energy functions to be used
      storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> > m_EnergyMap;

      //! set of SSMethods to use in evaluation
      storage::Set< sspred::Method> m_SSMethods;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief get the name of the object
      //! @return the name of the object
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PhiPsiWithSSPred();

      //! @brief constructor taking member variable parameters
      //! @param HISTOGRAM_FILENAME path to file where statistics and in consequence the energy potentials are read from
      //! @param SS_METHODS set of SSMethods to use in evaluation
      //! @param SCHEME scheme to be used in outputting
      PhiPsiWithSSPred
      (
        const storage::Set< sspred::Method> &SS_METHODS,
        const std::string &SCHEME = GetDefaultScheme(),
        const std::string &HISTOGRAM_FILENAME = PhiPsi::GetDefaultHistogramFilename()
      );

      //! @brief Clone function
      //! @return pointer to new PhiPsiWithSSPred
      PhiPsiWithSSPred *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the scheme
      //! @return string which is the scheme for outputting
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief gives the filename of the histogram where the statistics for phi psi scoring come from
      //! @return string which is the filename of the histogram where the statistics for phi psi scoring come from
      const std::string &GetHistogramFilename() const
      {
        return m_HistogramFileName;
      }

      //! @brief gives the energy map used to give the phi psi scores
      //! @return map of phi psi energies for each ss type and each aa type
      const storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> > &
      GetEnergyMap() const
      {
        return m_EnergyMap;
      }

      //! @brief gives the ss methods stored
      //! @return set of sspred methods which are used
      const storage::Set< sspred::Method> &GetSSMethods() const
      {
        return m_SSMethods;
      }

      //! @brief set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read from this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the score for a given SSE
      //! @param SSE SSE of interest
      //! @param MEMBRANE membrane object
      //! @return score calculated for the given SSE
      storage::Pair< double, size_t> operator()
      (
        const assemble::SSE &SSE,
        const biol::Membrane &MEMBRANE
      ) const;

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

    }; // class PhiPsiWithSSPred

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PHI_PSI_WITH_SSPRED_H_ 
