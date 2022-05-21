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

#ifndef BCL_SCORE_AA_PAIR_DISTANCE_FITTED_FUNCTION_H_
#define BCL_SCORE_AA_PAIR_DISTANCE_FITTED_FUNCTION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairDistanceFittedFunction
    //! @brief This is a Function derived class for storing the fitted scoring function to amino acid distance histograms
    //!
    //! @see @link example_score_aa_pair_distance_fitted_function.cpp @endlink
    //! @author woetzen
    //! @date 02.03.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairDistanceFittedFunction :
      public math::FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      storage::VectorND< 2, double> m_RepulsionStartXY;    //!< x and y of the repulsion where is starts
      storage::VectorND< 2, double> m_RepulsionEndXY;      //!< x and y of the repulsion where it ends
      storage::VectorND< 2, double> m_AttractionStartXY;   //!< x and y of the attraction where it starts
      storage::VectorND< 2, double> m_AttractionMinimumXY; //!< x and y of the attraction where it has its minimum
      storage::VectorND< 2, double> m_AttractionEndXY;     //!< x and y of the attraction where it ends

      //! default start values for the start of the repulsion term for the fitted aapairenergy distribution
      static const storage::VectorND< 2, double> s_DefaultRepulsionStartXY;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AAPairDistanceFittedFunction();

      //! @brief constructor from a given histogram
      //! @param AA_PAIR_DISTANCE_DISTRIBUTION Histogram that containes aa pair distance distribution
      AAPairDistanceFittedFunction
      (
        const math::Histogram &AA_PAIR_DISTANCE_DISTRIBUTION
      );

      //! @brief construct from a vector of bins and a given energy distribution
      //! @param BINNING binning to be used
      //! @param ENERGY_DISTRIBUTION energy distribution to be used
      AAPairDistanceFittedFunction
      (
        const linal::Vector< double> &BINNING,
        const linal::Vector< double> &ENERGY_DISTRIBUTION
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new AAPairDistanceFittedFunction copied from this one
      AAPairDistanceFittedFunction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to the end of the attraction, which is also the distance cutoff
      //! @return the distance above which the score will be 0
      double GetDistanceCutoff() const
      {
        return m_AttractionEndXY.First();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the energy according to the distance
      //! @param DISTANCE distance to be used
      //! @return the energy calculated for the given distance
      double operator()( const double &DISTANCE) const;

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
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculates the repulsion for the given distance
      //! @param DISTANCE distance to be used
      //! @return the repulsion for the given distance
      double Repulsion( const double &DISTANCE) const;

      //! @brief calculates the attraction for the given distance
      //! @param DISTANCE distance to be used
      //! @return the attraction for the given distance
      double Attraction( const double &DISTANCE) const;

    }; //class AAPairDistanceFittedFunction

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_PAIR_DISTANCE_FITTED_FUNCTION_H_
