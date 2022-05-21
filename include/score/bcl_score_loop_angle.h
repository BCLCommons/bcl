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

#ifndef BCL_SCORE_LOOP_ANGLE_H_
#define BCL_SCORE_LOOP_ANGLE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_energy_distribution.h"
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopAngle
    //! @brief This is a Function template derived class for scoring the loop angle between two secondary structure
    //! elements using a statistical derived potential, that bins the Euclidean distance between the ends of the SSEs
    //! (end of body) and the angle between the end of the SSEs and the center of gravity.
    //!
    //! @see @link example_score_loop_angle.cpp @endlink
    //! @author putnamdk, fischea, heinzes1
    //! @date Jan 17, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoopAngle :
      public ProteinModel // needs to be derived from score::ProteinModel, b/c it needs the protein's center of gravity
    {

    private:

    //////////
    // data //
    //////////

      std::string m_HistogramFileName; //! file from which statistics and in consequence energy potentials are read from
      size_t m_MaxSequenceDistance; //! maximum sequence distance between SSEs to be considered as consecutive
      std::string m_Scheme; //! scheme to be used in outputting schemes
      math::CubicSplineDamped m_EnergyFunctionShortLoops; //! energy function for short loops <= m_MaxSequenceDistance
      math::CubicSplineDamped m_EnergyFunctionLongLoops; //! energy function for long loops > m_MaxSequenceDistance

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where loop angle histograms and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default file where the loop angle table is stored
      //! @return default file where the loop angle table data is stored
      static const std::string &GetDefaultTableFilename();

      //! @brief returns the maximum sequence distance between SSEs to be considered as consecutive
      //! @return maximum sequence distance between SSEs to be considered as consecutive
      static size_t GetDefaultMaxmimumSequenceDistance();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief get the name of the object when used in a dynamic context
      //! @return  the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file and scheme
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param MAX_SEQ_DISTANCE maximum sequence distance between SSEs to be considered as consecutive
      //! @param SCHEME scheme to be used
      LoopAngle
      (
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const size_t &MAX_SEQ_DISTANCE = GetDefaultMaxmimumSequenceDistance(),
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new Loop copied from this one
      LoopAngle *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the histogram filename
      //! @return the histogram filename
      const std::string &GetHistogramFileName() const;

      //! @brief gets the maximum sequence distance for which the two SSEs are considered consecutive for this score
      //! @return the maximum sequence distance considered
      const size_t &GetMaxSequenceDistance() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get energy functions
      //! @return VectorND of short and long loop cubic spline
      storage::VectorND< 2, math::CubicSplineDamped> GetEnergyFunctions() const;

      //! @brief get score type
      //! @return score type
      ProteinModel::Type GetType() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the angle between the given SSEs respective to the given position
      //! @param SSE_FIRST first SSE of interest
      //! @param SSE_SECOND first SSE of interest
      //! @param CENTER_OF_MASS the center of mass of the protein the SSEs are in
      //! @return angle between the given SSEs respective to the given position
      double CalculateCosAngle
      (
        const assemble::SSE &SSE_FIRST,
        const assemble::SSE &SSE_SECOND,
        const linal::Vector3D &CENTER_OF_MASS
      ) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief the loop angle between the two given SSEs
      //! @param PROTEIN_MODEL the protein model to be scored
      //! @return score for the angle between consecutive SSEs
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
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write detailed scheme and values to OSTREAM
      //! @param MODEL protein models to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues( const assemble::ProteinModel &MODEL, std::ostream &OSTREAM) const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @return ERROR_STREAM stream with which to write errors
      bool ReadInitializeSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief read energy distribution for scoring pairs of AASequences
      void ReadEnergyVector();

      //! @brief helper function called by WriteDetailedSchemeAndValues and operator() so that the code remains in sync
      //! @param MODEL the protein model to be scored
      //! @param OSTREAM the output stream to write the detailed scheme to for this chain
      //! @param WRITE set to true to actually write to the output stream; otherwise, nothing will be written
      //! @return the final score
      double ScoreLoops( const assemble::ProteinModel &MODEL, std::ostream &OSTREAM, const bool &WRITE) const;

    }; // class LoopAngle

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_LOOP_ANGLE_H_
