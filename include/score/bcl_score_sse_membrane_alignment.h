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

#ifndef BCL_SCORE_SSE_MEMBRANE_ALIGNMENT_H_
#define BCL_SCORE_SSE_MEMBRANE_ALIGNMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEMembraneAlignment
    //! @brief This is a for scoring alignment of one SSE to the membrane normal
    //!
    //! @see @link example_score_sse_membrane_alignment.cpp @endlink
    //! @author woetzen
    //! @date 02.03.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEMembraneAlignment :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> >
    {

    private:

    //////////
    // data //
    //////////

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! energy map that stores the potentials for each SSType
      storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, math::CubicSplineDamped> > m_EnergyFunctions;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param SCHEME scheme to be used
      SSEMembraneAlignment
      (
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new SSEMembraneAlignment object copied from this one
      SSEMembraneAlignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access energy functions
      //! @return map for each sstype, that has a map for each environmenttype with a cubic spline
      const storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, math::CubicSplineDamped> > &
      GetEnergyFunctions() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the projection angle of the given SSE to the membrane plane
      //! @param SSE_GEOMETRY sse geometry derived class of interest
      //! @param AXIS membrane axis
      //! @param MEMBRANE membrane object
      //! @return calculated angle
      static double AngleToMembranePlane
      (
        const assemble::SSEGeometryInterface &SSE_GEOMETRY,
        const coord::Axis &AXIS,
        const biol::Membrane &MEMBRANE
      );

      //! @brief calculate weight for the given strand, that it is not turned in membrane
      //! @param STRAND sse geometry derived class of interest
      //! @param MEMBRANE membrane object
      //! @return calculated weight
      static double WeightXAxis
      (
        const assemble::SSEGeometryInterface &STRAND,
        const biol::Membrane &MEMBRANE
      );

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
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
      //! @param SSE SSE of interest
      //! @param MEMBRANE membrane object
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE,
        const biol::Membrane &MEMBRANE,
        std::ostream &OSTREAM
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief read energy distribution for scoring sse membrane alignment
      void
      ReadEnergyFunctions();

    }; //class SSEMembraneAlignment

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_SSE_MEMBRANE_ALIGNMENT_H_
