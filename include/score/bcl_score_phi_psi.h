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

#ifndef BCL_SCORE_PHI_PSI_H_
#define BCL_SCORE_PHI_PSI_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_energy_distribution.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PhiPsi
    //! @brief for scoring the agreement of phi psi values of residues an sse with expected probabilities
    //! @details For helix and strand SSEs, uses energy distributions specific to each of these two types of SSEs in
    //!          in order to determine the agreement of the argument SSE with the expected values of phi and psi within
    //!          either helix or strand. The type of amino acid is also taken into account
    //!          For coil, only the type of amino acid is taken into account. The ramachandran plots for each amino acid
    //!          type contain information for all SSE types.
    //!
    //!
    //! @see @link example_score_phi_psi.cpp @endlink
    //! @author rouvelgh, woetzen, karakam, alexanns, fischea
    //! @date Jun 30, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PhiPsi :
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

      //! map of energy functions to be used for soluble and membrane environments
      storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> > m_EnergyMapSoluble;
      storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> > m_EnergyMapMembrane;

      //! map of energy functions to be used for SSType independent phi psi scoring
      storage::Map< biol::AAType, math::BicubicSpline> m_AATypeEnergyMap;

      //! the types of sses whose phi psi angles will be scored independently of the sstype
      storage::Set< biol::SSType> m_SSTypes;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

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
      PhiPsi();

      //! @brief constructor from a specified histogram file
      //! @param SCHEME scheme to be used
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param SSTYPE_SET the types of sses whose phi psi angles will be scored independently of the sstype
      PhiPsi
      (
        const std::string &SCHEME,
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const storage::Set< biol::SSType> &SSTYPE_SET = storage::Set< biol::SSType>()
      );

      //! @brief Clone function
      //! @return pointer to new PhiPsi
      PhiPsi *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns filename of the histogram being used
      //! @return filename of the histogram being used
      const std::string &GetHistogramFilename() const
      {
        return m_HistogramFileName;
      }

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief access to the energy functions
      //! @return a map that has a 2D spline for phi and psi for each amino acid
      const storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> >
      &GetEnergyFunctions() const
      {
        return m_EnergyMapSoluble;
      }

      //! @brief access to the energy functions
      //! @return a map that has a 2D spline for phi and psi for each amino acid
      const storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> >
      &GetEnergyFunctionsMembrane() const
      {
        return m_EnergyMapMembrane;
      }

      //! @brief access to the energy functions for SSType independent phi psi scoring
      //! @return a map that has a 2D spline for phi and psi for each amino acid for SSType independent phi psi scoring
      const storage::Map< biol::AAType, math::BicubicSpline> &GetAATypeEnergyFunctions() const
      {
        return m_AATypeEnergyMap;
      }

      //! @brief gets the types of sses whose phi psi angles will be scores independently of the sstype
      //! @return the types of sses whose phi psi angles will be scores independently of the sstype
      const storage::Set< biol::SSType> &GetSSTypes() const
      {
        return m_SSTypes;
      }

      //! @brief set the ss types
      //! @param SS_TYPES ss types to be set
      void SetSSTypes( const storage::Set< biol::SSType> &SS_TYPES)
      {
        m_SSTypes = SS_TYPES;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief Operator that scores the phi_psi bending of various SSEs in a protein
      //! @param THIS_SSE SSE to be scored
      //! @param MEMBRANE membrane object
      //! @return pair of score and number of scored amino acids
      storage::Pair< double, size_t> operator()
      (
        const assemble::SSE &THIS_SSE,
        const biol::Membrane &MEMBRANE
      ) const;

      //! @brief scores an sse according to phi and psi without taking into account the SS Type
      //! @param THIS_SSE will be scored for agreement of phi and psi probabilities
      //! @return pair of double and size_t which is the score and number of scored amino acids, respectively
      storage::Pair< double, size_t> ScorePhiPsiSSTypeIndependent( const assemble::SSE &THIS_SSE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read from this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief read the energy distribution for scoring phi_psi angles
      void ReadEnergyVector();

      //! @brief scores the phi psi of a single residue according to a given energy map
      //! @param AATYPE_ENERGY the energy used to score the phi psi of the residue
      //! @param AA_BASE the residue of interest
      //! @param AA_PREVIOUS the residue previous in sequence to the residue of interest
      //! @param AA_NEXT the residue following in sequence to the residue of interest
      static double ScoreAAPhiPsi
      (
        const math::BicubicSpline &AATYPE_ENERGY,
        const biol::AABase &AA_BASE,
        const biol::AABase &AA_PREVIOUS,
        const biol::AABase &AA_NEXT
      );

      //! @brief scores the phi and psi angles of an entire aa sequence according to a given energy map
      //! @param AA_SEQUENCE the amino acid sequence which will be scored
      //! @param AATYPE_ENERGY_MAP the energy map used to score the amino acid sequence
      //! @return pair of double and size_t which are the score and number of residues scored, respectively
      static storage::Pair< double, size_t> ScoreSequencePhiPsi
      (
        const biol::AASequence &AA_SEQUENCE,
        const storage::Map< biol::AAType, math::BicubicSpline> &AATYPE_ENERGY_MAP_SOLUBLE,
        const storage::Map< biol::AAType, math::BicubicSpline> &AATYPE_ENERGY_MAP_MEMBRANE,
        const biol::Membrane &MEMBRANE
      );

    }; // class PhiPsi

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PHI_PSI_H_
