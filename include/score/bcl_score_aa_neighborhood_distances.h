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

#ifndef BCL_SCORE_AA_NEIGHBORHOOD_DISTANCES_H_
#define BCL_SCORE_AA_NEIGHBORHOOD_DISTANCES_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_neighborhood_interface.h"
#include "bcl_score_aa_pair_distance_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborhoodDistances
    //! @brief this class score the neighbors by adding all aa pair distance scores
    //!
    //! @see @link example_score_aa_neighborhood_distances.cpp @endlink
    //! @author woetzen
    //! @date 12.10.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborhoodDistances :
      public AANeighborhoodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! m_ScoreDistance as an implementation of AAPairDistanceInterface
      util::Implementation< AAPairDistanceInterface> m_ScoreDistance;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AANeighborhoodDistances();

      //! @brief constructor from a distance score
      //! @param SP_AA_PAIR_DISTANCE_SCORE aa pair distance score
      AANeighborhoodDistances
      (
        const AAPairDistanceInterface &SP_AA_PAIR_DISTANCE_SCORE
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new AANeighborhoodDistances object copied from this one
      AANeighborhoodDistances *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation for neighbors of that exposure score between amino acids in the same chain
      size_t GetMinimalSequenceSeparation() const;

      //! @brief access to the distance cutoff
      //! @return distance cutoff above which the neighbor does not have influence on the score anymore
      double GetDistanceCutoff() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate neighbor potential for a given amino acid and its AANeighborList
      //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
      //! @return neighbor potential for a given amino acid and its AANeighborList
      double operator()
      (
        const assemble::AANeighborList &AA_NEIGHBOR_LIST,
        const util::SiPtr< const biol::Membrane> &MEMBRANE
      ) const;

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
      //! @param INDENT number of indentations
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param AA_NEIGHBOR_LIST AANeighborList to be used for a single AABase
      //! @param OSTREAM the std::ostream to be written to
      //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::AANeighborList &AA_NEIGHBOR_LIST,
        const util::SiPtr< const biol::Membrane> &MEMBRANE,
        std::ostream &OSTREAM
      ) const;

      //! @brief return parameters for data members that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class AANeighborhoodDistances

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_NEIGHBORHOOD_DISTANCES_H_
