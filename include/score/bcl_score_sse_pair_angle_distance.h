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

#ifndef BCL_SCORE_SSE_PAIR_ANGLE_DISTANCE_H_
#define BCL_SCORE_SSE_PAIR_ANGLE_DISTANCE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_sse_pack_interface.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairAngleDistance
    //! @brief This is a wrapper class for using SSEPackInterface derived scoring classes within Pairs
    //!
    //! @see @link example_score_sse_pair_angle_distance.cpp @endlink
    //! @author karakam, woetzen
    //! @date 04/20/09
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPairAngleDistance :
      public math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double>
    {

    private:

    //////////
    // data //
    //////////

      //! ShPtr to scoring function that evaluates the the SSEGeometryPacking
      util::ShPtr< SSEPackInterface> m_SSEPackScoringFunction;

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
      SSEPairAngleDistance();

      //! @brief constructor from a SSEPackInterface
      //! @param SSE_PACK_INTERFACE reference to a SSEPackInterface derived class
      SSEPairAngleDistance( const SSEPackInterface &SSE_PACK_INTERFACE);

      //! @brief virtual copy constructor
      //! @return pointer to a new SSEPairAngleDistance copied from this one
      SSEPairAngleDistance *Clone() const;

      //! @brief destructor
      ~SSEPairAngleDistance();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the SSEGeometryPacking for given sse pair and return the score for this SSEGeometryPacking
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return the score for the calculated SSEGeometryPacking
      double operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
      {
        // if valid SSEs
        if( m_SSEPackScoringFunction->AreValidSSEs( SSE_A, SSE_B))
        {
          // calculate and return score
          return
          m_SSEPackScoringFunction->operator()
          (
            assemble::SSEGeometryPacking( SSE_A, SSE_B, m_SSEPackScoringFunction->GetMinimalInterfaceLength())
          );
        }
        // otherwise return 0
        else
        {
          return double( 0.0);
        }
      }

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
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::SSE &SSE_A,
        const assemble::SSE &SSE_B,
        std::ostream &OSTREAM
      ) const;

    }; //class SSEPairAngleDistance

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_SSE_PAIR_ANGLE_DISTANCE_H_
