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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_sse_geometry_packer_best_fragment_pairs.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackerBestFragmentPairs::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackerBestFragmentPairs())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPackerBestFragmentPairs::SSEGeometryPackerBestFragmentPairs() :
      m_Packer(),
      m_PackingCriteria(),
      m_PackingComparison()
    {
    }

    //! @brief constructor from a packer and a comparison function
    //! @param PACKER Packer to be used
    //! @param COMPARISION_FUNCTION function for comparing two SSE geometry packing objects
    SSEGeometryPackerBestFragmentPairs::SSEGeometryPackerBestFragmentPairs
    (
      const SSEGeometryPacker &PACKER,
      const math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> &COMPARISION_FUNCTION
    ) :
      m_Packer( PACKER),
      m_PackingCriteria(),
      m_PackingComparison( COMPARISION_FUNCTION.Clone())
    {
    }

    //! @brief constructor from a packer, a criteria function and a comparison function
    //! @param PACKER Packer to be used
    //! @param CRITERIA criteria to decide on which packings are to considered
    //! @param COMPARISION_FUNCTION function for comparing two SSE geometry packing objects
    SSEGeometryPackerBestFragmentPairs::SSEGeometryPackerBestFragmentPairs
    (
      const SSEGeometryPacker &PACKER,
      const math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> &CRITERIA,
      const math::BinaryFunctionInterface< SSEGeometryPacking, SSEGeometryPacking, bool> &COMPARISION_FUNCTION
    ) :
      m_Packer( PACKER),
      m_PackingCriteria( CRITERIA.Clone()),
      m_PackingComparison( COMPARISION_FUNCTION.Clone())
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEGeometryPackerBestFragmentPairs
    SSEGeometryPackerBestFragmentPairs *SSEGeometryPackerBestFragmentPairs::Clone() const
    {
      return new SSEGeometryPackerBestFragmentPairs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEGeometryPackerBestFragmentPairs::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the packing between all fragments of the given geometries and return it in a list
    //! @param SSE_GEOMETRY_A first SSEGeometryInterface derived class of interest
    //! @param SSE_GEOMETRY_B second SSEGeometryInterface derived class of interest
    //! @return list containing packing between all fragment pairings between given geometry pair
    storage::List< SSEGeometryPacking> SSEGeometryPackerBestFragmentPairs::operator()
    (
      const SSEGeometryInterface &SSE_GEOMETRY_A,
      const SSEGeometryInterface &SSE_GEOMETRY_B
    ) const
    {
      // static undefined packing
      static const SSEGeometryPacking s_undefined_packing;

      // if both SSEs are not helix or strand or they are the same geometry
      if
      (
        !SSE_GEOMETRY_A.GetType()->IsStructured() ||
        !SSE_GEOMETRY_B.GetType()->IsStructured() ||
        &SSE_GEOMETRY_A == &SSE_GEOMETRY_B
      )
      {
        // return empty list
        return storage::List< SSEGeometryPacking>();
      }

      // get the list of list of packings for the given two SSE Geometries
      storage::Vector< storage::List< SSEGeometryPacking> > all_packings_list
      (
        ( *m_Packer)->operator ()( SSE_GEOMETRY_A, SSE_GEOMETRY_B)
      );

      // initialize list to store the best fragments
      storage::List< SSEGeometryPacking> packing_list;

      // iterate over the first layers of lists
      for
      (
        storage::Vector< storage::List< SSEGeometryPacking> >::const_iterator list_itr( all_packings_list.Begin()),
          list_itr_end( all_packings_list.End());
        list_itr != list_itr_end; ++list_itr
      )
      {
        // create a pointer to store the best packing
        util::SiPtr< const SSEGeometryPacking> sp_best_packing( &s_undefined_packing);

        // iterate over the second layer of lists
        for
        (
          storage::List< SSEGeometryPacking>::const_iterator pack_itr( list_itr->Begin()),
          pack_itr_end( list_itr->End());
          pack_itr != pack_itr_end;
          ++pack_itr
        )
        {
          // boolean to see if the packing matches the criteria
          bool match_criteria( true);

          // if the criteria is set
          if( m_PackingCriteria.IsDefined())
          {
            match_criteria = m_PackingCriteria->operator()( *pack_itr);
          }

          // if this packing matches the criteria (if any) and
          // it is better than the best packing seen so far or this is the first one to be compared
          if
          (
            match_criteria &&
            ( !sp_best_packing->IsDefined() || m_PackingComparison->operator ()( *pack_itr, *sp_best_packing))
          )
          {
            // then update the best packing to this
            sp_best_packing = &( *pack_itr);
          }

        } // pack_itr
        // if the best packing is defined and has a valid interaction weight
        if( sp_best_packing->IsDefined())
        {
          // insert this list into the main list
          packing_list.PushBack( *sp_best_packing);
        }
      } // list_itr

      // end
      return packing_list;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackerBestFragmentPairs::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Packer, ISTREAM);
      io::Serialize::Read( m_PackingCriteria, ISTREAM);
      io::Serialize::Read( m_PackingComparison, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackerBestFragmentPairs::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Packer, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PackingCriteria, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PackingComparison, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
