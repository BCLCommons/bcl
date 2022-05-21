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

#ifndef BCL_ASSEMBLE_PICK_SSE_FURTHEST_EUCLIDEAN_CENTER_H_
#define BCL_ASSEMBLE_PICK_SSE_FURTHEST_EUCLIDEAN_CENTER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_pick_sse_furthest_euclidean.h"
#include "bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickSSEFurthestEuclideanCenter
    //! @brief class is used for picking the single furthest item from a point of reference.
    //! @details This class picks from a list of SSEs, the one that is furthest from the center of the given criteria
    //!
    //! @see @link example_assemble_pick_sse_furthest_euclidean_center.cpp @endlink
    //! @author karakam
    //! @date 04/23/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_CriteriaType>
    class PickSSEFurthestEuclideanCenter :
      public find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, t_CriteriaType>
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PickSSEFurthestEuclideanCenter()
      {
      }

      //! virtual copy constructor
      PickSSEFurthestEuclideanCenter< t_CriteriaType> *Clone() const
      {
        return new PickSSEFurthestEuclideanCenter< t_CriteriaType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! Picks the assemble::SSE object which is furthest away from "m_ReferencePoint"
      //! @param SSE_LIST list of SSEs
      //! @param CRITERIA criteria
      //! @return returns SiPtr to the assemble::SSE object which is furthest away from "m_ReferencePoint"
      util::SiPtr< const SSE>
      Pick( const util::SiPtrList< const SSE> &SSE_LIST, const t_CriteriaType &CRITERIA) const
      {
        // create PickSSEFurthestEuclidean and pick with the ProteinModelCenter
        return PickSSEFurthestEuclidean().Pick( SSE_LIST, CRITERIA.GetCenter());
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT indentation
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    }; // template class PickSSEFurthestEuclideanCenter

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PICK_SSE_FURTHEST_EUCLIDEAN_CENTER_H_
