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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_BEST_FRAGMENT_PAIR_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_BEST_FRAGMENT_PAIR_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_packers.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackerBestFragmentPair
    //! @brief class returns the best packing for a given geometry pair
    //! @details This class uses a member SSEGeometryPacker to get list of list of SSEGeometryPackings between all the
    //! sub-geometries for the given pair of SSEGeometryInterface derived classes. It uses the members criteria
    //! and the comparison operator to decide on the best packing and returns it.
    //! For GEOMETRY_A ( composed of M sub-geometries) and GEOMETRY_B (composed of N sub-geometry)
    //! only one packing will be returned
    //!
    //! @see @link example_assemble_sse_geometry_packer_best_fragment_pair.cpp @endlink
    //! @author karakam
    //! @date Jul 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackerBestFragmentPair :
      public SSEGeometryPackerInterface< SSEGeometryPacking>
    {

    private:

    //////////
    // data //
    //////////

      //! function for calculating the actual packings between fragments
      SSEGeometryPacker m_Packer;

      //! criteria to decide on which packings are to considered
      util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > m_PackingCriteria;

      //! function for comparing two SSE geometry packing objects
      util::ShPtr< math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> > m_PackingComparison;

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
      SSEGeometryPackerBestFragmentPair();

      //! @brief constructor from a packer and a comparison function
      //! @param PACKER SSEGeometryPacker to be used
      //! @param COMPARISION_FUNCTION function for comparing two SSE geometry packing objects
      SSEGeometryPackerBestFragmentPair
      (
        const SSEGeometryPacker &PACKER,
        const math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> &COMPARISION_FUNCTION
      );

      //! @brief constructor from a packer, a criteria function and a comparison function
      //! @param PACKER SSEGeometryPacker to be used
      //! @param CRITERIA criteria to decide on which packings are to considered
      //! @param COMPARISION_FUNCTION function for comparing two SSE geometry packing objects
      SSEGeometryPackerBestFragmentPair
      (
        const SSEGeometryPacker &PACKER,
        const math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> &CRITERIA,
        const math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> &COMPARISION_FUNCTION
      );

      //! @brief Clone function
      //! @return pointer to new SSEGeometryPackerBestFragmentPair
      SSEGeometryPackerBestFragmentPair *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the minimal interface length used
      //! @return minimal interface length used
      double GetMinimalInterfaceLength() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the packing between all fragments of the given geometries and return the best one
      //! @param SSE_GEOMETRY_A first SSEGeometryInterface derived class of interest
      //! @param SSE_GEOMETRY_B second SSEGeometryInterface derived class of interest
      //! @return the best packign between fragments of given geometry pair according to comparison function
      SSEGeometryPacking operator()
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B
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

    }; // class SSEGeometryPackerBestFragmentPair

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_BEST_FRAGMENT_PAIR_H_
