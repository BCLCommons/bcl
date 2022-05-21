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

#ifndef BCL_SCORE_LOOP_H_
#define BCL_SCORE_LOOP_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Loop
    //! @brief This is a Function derived template class for scoring the loop length between two secondary structure elements
    //! using a statistical derived potential, that bin the sequence distance (different is seq id of end of first sse
    //! and begin of second sse) and the Euclidean distance between the ends of the sse (end of body).
    //!
    //! @see @link example_score_loop.cpp @endlink
    //! @author woetzen
    //! @date 10.06.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API Loop :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! scheme to be used in outputting schemes
      std::string m_Scheme;

      //! max number of residues in loop
      size_t m_MaxLoopLength;

      //! energy function to be used binned by nr residues up to max loop length - everything larger will be scored with
      storage::Vector< math::CubicSplineDamped> m_EnergyFunctions;

      //! energy function for longer loops than m_MaxLoopLength
      math::CubicSplineDamped m_EnergyFunctionLongLoops;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns default file where loop distance histograms and in consequence energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default file where the loop distance table is stored
      //! @return default file where the loop distance table data is stored
      static const std::string &GetDefaultTableFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief returns default maximum loop length in residues to be considered
      //! @return default maximum loop length in residues to be considered
      static size_t GetDefaultMaxLoopLength();

      //! @brief get the name of the object
      //! @return the name of the object
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param SCHEME scheme to be used
      //! @param MAX_LOOP_LENGTH maximum number of residues in a loop to be considered
      Loop
      (
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const std::string &SCHEME = GetDefaultScheme(),
        const size_t MAX_LOOP_LENGTH = GetDefaultMaxLoopLength()
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new Loop copied from this one
      Loop *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief get maximal loop length, above which there is only one energy distribution for Euclidean distance over log( seq distance + 2)
      //! @return max loop length
      size_t GetMaxLoopLength() const
      {
        return m_MaxLoopLength;
      }

      //! @brief return reference to the energy function
      //! @return const ref to a CubicSpline where x is the seq distance angle and y the Euclidean distance
      const storage::Vector< math::CubicSplineDamped> &GetEnergyFunctions() const
      {
        return m_EnergyFunctions;
      }

      //! @brief get energy function above max loop length
      //! @return cubic spline that evaluates Euclidean distance over log( seq distance + 2)
      const math::CubicSplineDamped &GetEnergyFunctionAboveMaxLoopLength() const
      {
        return m_EnergyFunctionLongLoops;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief score the loop length and distance
      //! @param SEQ_EUC_DISTANCE pair of sequence and Euclidean distance
      //! @return score for that combination of sequence and Euclidean distance
      double Score( const storage::Pair< size_t, double> &SEQ_EUC_DISTANCE) const;

      //! @brief calculate Sequence distance between aas and Euclidean distance between ends of bodies of two SSEs
      //! @brief SSE_A first SSE of interest
      //! @brief SSE_B first SSE of interest
      //! @return a pair which first member is the length of the loop and aminoacids and the second is the Euclidean distance
      static
      storage::Pair< size_t, double>
      SequenceAndEuclideanDistance
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B
      );

      //! @brief normalizes the distance by the sequence length
      //! @param SEQ_EUC_DISTANCE pair of sequence and Euclidean distance
      //! @return normalized distance
      static double NormalizeDistance
      (
        const storage::Pair< size_t, double> &SEQ_EUC_DISTANCE
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief the loop length and distance between the two given sses
      //! @brief SSE_A first SSE of interest
      //! @brief SSE_B first SSE of interest
      //! @return potential for loop
      double operator()
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    public:

      //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
      //! @param SSE_A first SSE of interest
      //! @param SSE_B first SSE of interest
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B,
        std::ostream &OSTREAM
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculate the maximum loop euclidian distance from the histogram files for each loop sequence length within a certain percentile
      //! @param HISTOGRAM_FILE filename of histogram to be used
      //! @param MAX_NR_LOOP_RESIDUES maximum number of loop residues
      //! @param FRACTION fraction of counts that needs to be below this threshold
      static storage::Vector< double> CalculateMaximumObservedDistances
      (
        const std::string &HISTOGRAM_FILENAME,
        const size_t MAX_NR_LOOP_RESIDUES,
        const double FRACTION
      );

      //! @brief set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read from this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    private:

      //! @brief read energy distribution for scoring pairs of AASequences
      void ReadEnergyVector();

    }; //class Loop

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_LOOP_H_
