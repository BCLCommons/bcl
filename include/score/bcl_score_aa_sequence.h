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

#ifndef BCL_SCORE_AA_SEQUENCE_H_
#define BCL_SCORE_AA_SEQUENCE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_aa_pair_distance_interface.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASequence
    //! @brief FunctionInterface derived class for scoring amino acids within a single SSE
    //! @details This is a FunctionInterface derived class for scoring a single SSE using a pairwise
    //! amino acid interaction potential.
    //!
    //! @see @link example_score_aa_sequence.cpp @endlink
    //! @author karakam, alexanns
    //! @date 16.01.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASequence :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> >
    {

    private:

    //////////
    // data //
    //////////

      //! Function to score amino acid pair interaction
      util::ShPtr< AAPairDistanceInterface> m_ScoreAAPair;

      //! normalize by number of scored aa pairs (score != 0)
      bool m_Normalize;

      //! function to generate the neighbor list
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::SSE, assemble::AANeighborListContainer> > m_AANeighborListContainerGenerator;

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
      AASequence();

      //! @brief construct from ShPtr to amino acid pair potential
      //! @param SP_AA_PAIR_POTENTIAL ShPtr to the amino acid potential to be used
      //! @param NORMALIZE if true, score will be normalize by number of non zero pair scores
      AASequence
      (
        const util::ShPtr< AAPairDistanceInterface> &SP_AA_PAIR_POTETNIAL,
        const bool NORMALIZE
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new AASequence copied from this one
      AASequence *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief score aa distances between all pairs of amino acids within a single SSE
      //! @param THIS_SSE SSE of interest
      //! @param MEMBRANE membrane object
      //! @return pair of overall interaction potential and number of scored entities
      storage::Pair< double, size_t> operator()
      (
        const assemble::SSE &THIS_SSE,
        const biol::Membrane &MEMBRANE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! read from std::istream
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
      //! @param THIS_SSE SSE of interest
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::SSE &THIS_SSE,
        std::ostream &OSTREAM
      ) const;

    }; // class AASequence

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_SEQUENCE_H_
