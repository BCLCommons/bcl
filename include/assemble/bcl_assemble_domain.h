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

#ifndef BCL_ASSEMBLE_DOMAIN_H_
#define BCL_ASSEMBLE_DOMAIN_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_domain_interface.h"
#include "bcl_assemble_sse.h"
#include "bcl_assemble_sse_compare.h"
#include "bcl_assemble_topology.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Domain
    //! @brief contains Set of SSEs and provides functionality to manipulate them
    //! @details Set of SSEs are stored by this class in a sorted fashion where no overlapping SSEs can be found. It
    //! provides mutate functions as well as insertion, removal, replacement for SSEs. This class forms the SSE part of
    //! a chain. It has no direct information about the full sequence, it only knows the sequence parts represented by
    //! the SSEs stored in the set. It can also be used store SSEs from multiple chains when used separately.
    //!
    //! @see @link example_assemble_domain.cpp @endlink
    //! @author woetzen, karakam
    //! @date Nov 5, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Domain :
      public DomainInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Set of ShPtr to SSEs that form domain
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> m_Data;

      //! Topology
      util::ShPtr< Topology> m_Topology;

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
      Domain();

      //! @brief construct from a Set of ShPtr to SSEs
      //! @param SSE_SET Set of ShPtr to SSEs
      Domain( const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &SSE_SET);

      //! @brief construct from a Set of ShPtr to SSEs and a corresponding topology
      //! @param SSE_SET Set of ShPtr to SSEs
      //! @param SP_TOPOLOGY ShPtr to corresponding Topology
      Domain
      (
        const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &SSE_SET,
        const util::ShPtr< Topology> &SP_TOPOLOGY
      );

      //! @brief construct from util::ShPtrVector of SSE
      //! @param SSE_VECTOR ShPtrVector of SSEs
      Domain( const util::ShPtrVector< SSE> &SSE_VECTOR);

      //! @brief construct from util::ShPtrVector of SSE and a corresponding topology
      //! @param SSE_VECTOR ShPtrVector of SSEs
      //! @param SP_TOPOLOGY ShPtr to corresponding Topology
      Domain
      (
        const util::ShPtrVector< SSE> &SSE_VECTOR,
        const util::ShPtr< Topology> &SP_TOPOLOGY
      );

      //! @brief copy constructor for domain
      //! @param DOMAIN_RHS Domain to be copied
      Domain( const Domain &DOMAIN_RHS);

      //! @brief virtual copy constructor
      //! @return pointer to a new Domain instance copied from this one
      Domain *Clone() const;

      //! @brief virtual hard copy constructor
      //! @return pointer to a new Domain instance hard-copied from this one
      Domain *HardCopy() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return all chain ids within that domain
      //! @return set of all chain ids
      storage::Set< char> GetChainIds() const;

      //! @brief return total number of sses
      //! @return total number of sses
      size_t GetNumberSSEs() const
      {
        return m_Data.GetSize();
      }

      //! @brief return number of SSE of specified SSTYPE
      //! @param SS_TYPE specific SSTYPE
      //! @return number of SSE of specified SSTYPE
      size_t GetNumberSSE( const biol::SSType &SS_TYPE) const;

      //! @brief returns the number of amino acids in the chain
      //! @return the number of amino acids in the chain
      size_t GetNumberAAs() const;

      //! @brief return const m_Data
      //! @return const m_Data
      const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &GetData() const
      {
        return m_Data;
      }

      //! @brief returns all SSEs in domain in a util::SiPtrVector
      //! @return all SSEs in domain in a util::SiPtrVector
      util::SiPtrVector< const SSE> GetSSEs() const;

      //! @brief returns all SSEs in domain of given SSTYPE in a util::SiPtrVector
      //! @param SS_TYPE specific SSTYPE
      //! @return all SSEs in domain of given SSTYPE in a util::SiPtrVector
      util::SiPtrVector< const SSE> GetSSEs( const biol::SSType &SS_TYPE) const;

      //! @brief set m_Data to provided SSE_SET
      //! @param SSE_SET set of SSEs
      void SetData( const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &SSE_SET)
      {
        m_Data = SSE_SET;
      }

      //! @brief return topology associated with this domain
      //! @return topology associated with this domain
      const util::ShPtr< Topology> &GetTopology() const
      {
        return m_Topology;
      }

      //! @brief set topology to the given one
      //! @param SP_TOPOLOGY new topology to be used
      void SetTopology( const util::ShPtr< Topology> &SP_TOPOLOGY)
      {
        m_Topology = SP_TOPOLOGY;
      }

      //! @brief concatenates sequences of all sses fills in the gaps with unknown AAs and returns this AASequence
      //! @return AASequence created by concatenating sequences of all sses and filling in the gaps with unknown AAs
      util::ShPtr< biol::AASequence> CreateSequenceFromSSEs() const;

      //! @brief returns the geometric center of the object
      //! @return the geometric center of the object
      linal::Vector3D GetCenter() const;

      //! @brief return the orientation of the object
      //! @return orientation
      linal::Vector3D GetAxis( const coord::Axis &AXIS) const;

      //! @brief return the orientation and Position as TransformationMatrix3D
      //! @return TransformationMatrix3D that defines orientation and position
      const math::TransformationMatrix3D GetOrientation() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief translates the coordinates of all SSEs by the supplied translation vector
      //! @param TRANSLATION_VECTOR_3D Translation vector to be applied
      void Translate( const linal::Vector3D &TRANSLATION_VECTOR_3D);

      //! @brief transforms the coordinates of all SSEs according to given transformation matrix
      //! @param TRANSFORMATION_MATRIX_3D transformation matrix to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief rotate the SSE by a given rotation matrix
      //! @param ROTATION_MATRIX_3D rotation matrix to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

      //! @brief pushback a util::ShPtr< SSE> to m_Data
      //! @param SSELEMENT ShPtr to SSE to be inserted
      //! @return whether insertion succeeded
      bool Insert( const util::ShPtr< SSE> &SSELEMENT);

      //! @brief pushback a vector SSEs
      //! @param SSE_VECTOR Vector of SSEs to be added
      //! @return whether insertion succeeded
      bool Insert( const util::ShPtrVector< SSE> &SSE_VECTOR);

      //! @brief insert the SSEs from the given domain
      //! @param NEW_DOMAIN Domain from which SSEs should be inserted
      //! @return whether insertion was successful
      bool Insert( const Domain &NEW_DOMAIN);

      //! @brief replace the given SP_SSE with already existing one
      //! @param SP_SSE ShPtr pointing to the SSE to be replaced
      bool Replace( const util::ShPtr< SSE> &SP_SSE);

      //! @brief replace all SSEs that overlap with SP_SSE with SP_SSE
      //! @param SP_SSE ShPtr to SSE to be inserted
      //! @return whether replacement succeeded
      bool ReplaceWithOverlapping( const util::ShPtr< SSE> &SP_SSE);

      //! @brief remove given SSELEMENT from the domain
      //! @param SSELEMENT SSE to be removed
      bool Remove( const SSE &SSELEMENT);

      //! @brief sets positions of all SSEs to ideal conformation w/wo superimposing with prior coordinates
      //! @param KEEP_POSITION flag to indicate whether to original body information of SSEs should be reserved
      void SetToIdealConformation( const bool KEEP_POSITION = true);

      //! @brief chop all SSE elements of that model in pieces of the sizes defined by MIN_SSE_LENGTHS
      //! @param MIN_SSE_LENGTHS VectorND of sizes that defined min size for each SSType
      void ChopSSEs( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS);

      //! @brief find and to return the ShPtr for the given SSE
      //! @param SSE_TO_SEARCH SSE of interest
      //! @return ShPtr to corresponding SSE, otherwise an empty ShPtr
      const util::ShPtr< SSE> &FindSSE( const SSE &SSE_TO_SEARCH) const;

      //! @brief checks if domain already contains THIS_SSE
      //! @param THIS_SSE SSE of interest
      //! @return whether domain already contains THIS_SSE
      bool DoesContain( const SSE &THIS_SSE) const;

      //! @brief returns true if the domain has no SSEs in it
      //! @return true if the domain has no SSEs in it
      bool IsEmpty() const
      {
        return m_Data.IsEmpty();
      }

      //! @brief checks if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
      //! @param THIS_SSE SSE to be searched for
      //! @return if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
      bool DoesContainOverlapping( const SSE &THIS_SSE) const;

      //! @brief join following ( progressing sequence id) SSEs of given SS_TYPE into one SSE
      //! @param SS_TYPE SSType of interest
      void Join( const biol::SSType &SS_TYPE);

      //! @brief filters the current chain by given minimum SSE sizes
      //! @param MIN_SSE_SIZES minimum SSE sizes to filter the chain by
      void FilterByMinSSESizes( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES);

    ///////////////
    // operators //
    ///////////////

      //! @brief equal operator
      //! @param DOMAIN_RHS Domain to be assigned to
      //! @return this domain after being assigned to DOMAIN_RHS
      Domain &operator =( const Domain &DOMAIN_RHS);

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read Domain from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write Domain to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class Domain

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_DOMAIN_H_
