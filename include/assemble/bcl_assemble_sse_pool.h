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

#ifndef BCL_ASSEMBLE_SSE_POOL_H_
#define BCL_ASSEMBLE_SSE_POOL_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "command/bcl_command.fwd.hh"
#include "pdb/bcl_pdb.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_domain_interface.h"
#include "bcl_assemble_sse.h"
#include "bcl_assemble_sse_compare.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPool
    //! @brief This class stores possible sses derived from predictions in an organized structure
    //! @details This class provides storage for different SSEs that can be derived from an SSE. It's produced by any
    //! SSEFactoryInterface derived class and the stored SSEs are allowed to overlap with each other. It also provides
    //! convenience functions to read and write in PDB-like SSE definition format
    //!
    //! @see @link example_assemble_sse_pool.cpp @endlink
    //! @author karakam, woetzen, linders
    //! @date 12.02.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPool :
      public DomainInterface
    {

    private:

    //////////
    // data //
    //////////

      //! set that contains all SSEs of the pool
      storage::Set< util::ShPtr< SSE>, SSELessThan> m_Data;

    public:

      //! @brief const_iterator on data
      typedef storage::Set< util::ShPtr< SSE>, SSELessThan>::const_iterator const_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief return command line flag for providing a pool file
      //! @return command line flag for providing a pool file
      static util::ShPtr< command::FlagInterface> &GetFlagPoolRead();

      //! @brief return command line flag for setting minimum helix and strand lengths to be considered from the pool
      //! @return command line flag setting minimum helix and strand lengths to be considered from the pool
      static util::ShPtr< command::FlagInterface> &GetFlagMinSSELengths();

      //! @brief return command line flag for setting minimum helix and strand lengths to be considered from the pool
      //! @return command line flag setting minimum helix and strand lengths to be considered from the pool
      static storage::Map< biol::SSType, size_t> GetCommandLineMinSSELengths();

      //! @brief return command line flag for specifying a pool prefix - to be used when different pools postfixes are
      //!        specified within a stages file
      //! @return command line flag setting the prefix of pool files that will be used
      static util::ShPtr< command::FlagInterface> &GetFlagPoolPrefix();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEPool();

      //! @brief constructor from a vector of sses
      //! @param SSE_VECTOR vector of sses
      //! @param IGNORE_UNSTRUCTURED ignore unstructured sses
      //! @param IDEALIZE whether to idealize the SSEs (or simply move them to the origin)
      SSEPool
      (
        const util::SiPtrVector< const SSE> &SSE_VECTOR,
        const bool IGNORE_UNSTRUCTURED = true,
        const bool IDEALIZE = true
      );

      //! @brief constructor from a list of sses
      //! @param SSE_LIST list of sses
      //! @param IGNORE_UNSTRUCTURED ignore unstructured sses
      //! @param IDEALIZE whether to idealize the SSEs (or simply move them to the origin)
      SSEPool
      (
        const util::SiPtrList< const SSE> &SSE_LIST,
        const bool IGNORE_UNSTRUCTURED = true,
        const bool IDEALIZE = true
      );

      //! @brief virtual copy constructor
      SSEPool *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns SiPtrVector of all SSEs
      //! @return SiPtrVector of all SSEs
      util::SiPtrVector< const SSE> GetSSEs() const;

      //! @brief returns SiPtrVector of all SSEs in pool of given SS_TYPE
      //! @param SS_TYPE SSType of interest
      //! @return SiPtrVector of all SSEs of given SS_TYPE
      util::SiPtrVector< const SSE> GetSSEs( const biol::SSType &SS_TYPE) const;

      //! @brief Get the total number of potentially structured AAs
      size_t GetNumberPotentiallyStructuredAAs() const;

      //! @brief Get the total number of potentially structured AAs
      size_t GetNumberAAs() const
      {
        return GetNumberPotentiallyStructuredAAs();
      }

      //! @brief number of sses in pool
      //! @return size of pool
      size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      //! @brief iterator on begin
      //! @return const_iterator pointing to begin
      const_iterator Begin() const
      {
        return m_Data.Begin();
      }

      //! @brief iterator on end
      //! @return const_iterator pointing to end
      const_iterator End() const
      {
        return m_Data.End();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief insert sses into pool
      //! @param ITR_BEGIN iterator on first element of sses to be inserted
      //! @param ITR_END iterator to end of sses to be insered
      //! the iterators should be ShPtr< SSE>
      template< typename t_IteratorType>
      void InsertElements( const t_IteratorType &ITR_BEGIN, const t_IteratorType &ITR_END)
      {
        m_Data.InsertElements( ITR_BEGIN, ITR_END);
      }

      //! @brief insert sse pool's sses into pool
      //! @param SSE_POOL pool of sses, which sses should be inserted
      void InsertElements( const SSEPool &SSE_POOL)
      {
        m_Data.InsertElements( SSE_POOL.m_Data);
      }

      //! @brief insert sse into pool
      //! @param SP_SSE sse to be inserted
      void Insert( const util::ShPtr< SSE> &SP_SSE)
      {
        m_Data.Insert( SP_SSE);
      }

      //! @brief returns whether the pool has overlapping SSEs
      //! @return whether the pool has overlapping SSEs
      bool IsOverlapping() const;

      //! @brief removes the identical SSEs with the provided PROTEIN_MODEL from pool
      //! @param PROTEIN_MODEL protein model
      //! @return SiPtrList of SSEs in the pool excluding ones that are identical to SSEs in the provided PROTEIN_MODEL
      util::SiPtrList< const SSE> GetNonIdenticalSSEs
      (
        const ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief removes the overlapping SSEs with the provided PROTEIN_MODEL from pool
      //! @param PROTEIN_MODEL protein model
      //! @return SiPtrList of SSEs where there are no overlapping SSEs with the provided PROTEIN_MODEL
      util::SiPtrList< const SSE> GetNonOverlappingSSEs( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief return the SSEs of the given type in the pool which do not overlap with SSEs in the given model
      //! @param PROTEIN_MODEL protein model to compare to
      //! @param SSE_TYPE which SSE type to consider
      //! @return SSEs of the given type in the pool which do not overlap with SSEs in the given model
      util::SiPtrList< const SSE> GetNonOverlappingSSEs
      (
        const ProteinModel &PROTEIN_MODEL,
        const biol::SSType &SSE_TYPE
      ) const;

      //! @brief returns overlapping SSEs in the pool with the SSEs from provided SSE_DOMAIN
      //! @param SSE_DOMAIN domain
      //! @param EXCLUDE_IDENTICAL boolean to decide whether SSEs from domain should be included in the return
      //! @param EXCLUDE_DIFFERENT_SSTYPE whether to exclude overlapping SSEs that differ in SSType
      //! @return SiPtrList of SSEs which overlap with SSEs with the provided SSE_DOMAIN
      util::SiPtrList< const SSE> GetOverlappingSSEs
      (
        const DomainInterface &SSE_DOMAIN,
        const bool EXCLUDE_IDENTICAL = true,
        const bool EXCLUDE_DIFFERENT_SSTYPE = false
      ) const;

      //! @brief returns overlapping SSEs in the pool with the provided TARGET_SSE
      //! @param TARGET_SSE SSE of interest
      //! @param EXCLUDE_IDENTICAL boolean to decide whether SSEs from ProteinModel should be included in the return
      //! @param EXCLUDE_DIFFERENT_SSTYPE whether to exclude overlapping SSEs that differ in SSType
      //! @return SiPtrList of SSEs which overlap with TARGET_SSE
      util::SiPtrList< const SSE> GetOverlappingSSEs
      (
        const SSE &TARGET_SSE,
        const bool EXCLUDE_IDENTICAL = true,
        const bool EXCLUDE_DIFFERENT_SSTYPE = false
      ) const;

      //! @brief converts the data in the set from SSELessThan to SSELessThanNoOverlap
      //! @return SSELessThanNoOverlap set
      storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> GetRandomNonOverlappingSet() const;

      //! @brief calculates the number of helices and strands in the pool
      //!        this will equal to actual counts if it's non overlapping pool
      //!        if overlapping, it will be average numbers that lead to a complete model with non-overlapping SSEs
      //! @return average number of helices and strands expected in a complete model from this pool
      storage::Pair< double, double> CalculateAverageHelixStrandCounts() const;

      //! @brief calculate an estimate of the helix to strand ratio and return it
      //! @return helix to strand ratio ( 0:1 if only strands, 1:0 if only helices)
      storage::Pair< double, double> CalculateHelixToStrandRatio() const;

      //! @brief number of sses above given length
      //! @param MIN_SIZE_HELIX
      //! @return nr strands size > MIN_SIZE_STRAND + nr helices size > MIN_HELIX_SIZE
      double CountLongNonOverlappingSSEs( const size_t MIN_SIZE_STRAND, const size_t MIN_SIZE_HELIX) const;

      //! @brief prune the pool from sses that are not in the given map or do not have the minimum length
      //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
      void Prune( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP);

      //! @brief chops the SSEs in of the type in the map into chunks of the minimal size given in the map
      //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
      //! @return true, if any sses were joined, false otherwise
      void ChopSSEs( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP);

      //! @brief join adjacent sses if the adjacent sse is shorter than the given sse min size and has the same type, starting from the outside of a sequence of adjacent sses
      //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
      //! @return true if at least one join was performed - one call might not perform all possible joins and it can be necessary to call that function multiple timess
      bool Join( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP);

      //! @brief separate adjacent sses of identical type, by 2 * nr residues, but only if the resulting sses still have the desired size
      //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
      //! @param NR_RESIDUES number of residues
      //! @return true if at least one separation was performed
      bool Separate( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP, const size_t NR_RESIDUES);

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION Translation to be applied
      void Translate( const linal::Vector3D &TRANSLATION)
      {
      }

      //! @brief transform the object by a given TransformationMatrix3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
      {
      }

      //! @brief rotate the object by a given RotationMatrix3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
      {
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read SSEPool from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write SSEPool to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief read SSEPool from std::istream from a specially formatted pool file
      //! @param ISTREAM input stream
      //! @param PROTEIN_MODEL ProteinModel which SSEs belong to
      //! @param MIN_HELIX_LENGTH Minimal length for a HELIX to be added to the pool
      //! @param MIN_STRAND_LENGTH Minimal length for a STRAND to be added to the pool
      //! @return istream which was read from
      std::istream &ReadSSEPool
      (
        std::istream &ISTREAM,
        const ProteinModel &PROTEIN_MODEL,
        const size_t MIN_HELIX_LENGTH = 7,
        const size_t MIN_STRAND_LENGTH = 5
      );

      //! @brief read SSEPool from std::istream from a specially formatted pool file
      //! @param ISTREAM input stream
      //! @param SEQUENCE aa sequence to which SSEs belong to
      //! @param MIN_HELIX_LENGTH Minimal length for a HELIX to be added to the pool
      //! @param MIN_STRAND_LENGTH Minimal length for a STRAND to be added to the pool
      //! @return istream which was read from
      std::istream &ReadSSEPool
      (
        std::istream &ISTREAM,
        const biol::AASequence &SEQUENCE,
        const size_t MIN_HELIX_LENGTH = 7,
        const size_t MIN_STRAND_LENGTH = 5
      );

      //! @brief ReadSSEPoolInformation reads the pool file and gives the information it contains
      //!        This function does not create an SSE pool.
      //! @param POOL_FILENAME is the name and path of the file which contains the SSE pool
      //! @return returns a storage::Vector of the pdb::Lines which are contained in the SSE pool
      static storage::Vector< pdb::Line> ReadSSEPoolInformation
      (
        const std::string &POOL_FILENAME
      );

      //! @brief GetChainsRepresented obtains a storage::Set of all the chain ids in a pool file
      //!        This function does not create an SSE pool.
      //! @param POOL_FILENAME is the name and path of the file which contains the SSE pool
      //! @return returns a storage::Set< char> which has all the chain ids contained in the SSE pool
      static storage::Set< char> GetChainsRepresented
      (
        const std::string &POOL_FILENAME
      );

    public:

      //! @brief write SSEPool to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &WriteSSEPool( std::ostream &OSTREAM) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief finds SSEs in pool that agree best with SSEs in PROTEIN_MODEL
      //! @param PROTEIN_MODEL protein model whose SSEs are to be checked against the pool
      //! @return ShPtrVector to SSEs in pool that agree best with SSEs in PROTEIN_MODEL
      util::SiPtrVector< const SSE> FindBestMatchFromPool( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief finds SSEs in pool that agree best with SSE that is passed to function
      //! @param TARGET_SSE SSE to be checked against the pool
      //! @param MATCH_THRESHOLD the multiplier for setting the minimum threshold to be considered a match. This value
      //!                        is multiplied by the length of the TARGET_SSE to get the minimum threshold value
      //! @return ShPtr to SSE in pool that agrees best with SSE that is passed to function
      storage::Pair< util::SiPtr< const SSE>, double>
      FindBestMatchFromPool( const SSE &TARGET_SSE, const double MATCH_THRESHOLD = double( 2)) const;

      //! @brief calculates overlap measure between TARGET_SSE and SSE_FROM_POOL
      //! @param TARGET_SSE SSE to be checked against SSE from pool
      //! @param SSE_FROM_POOL SSE to be checked against TARGET_SSE
      //! @return int overlap measure
      int CalculateOverlapMeasure( const SSE &TARGET_SSE, const SSE &SSE_FROM_POOL) const;

      //! @brief function to return statistics table header names
      //! @return vector of table header names
      static const storage::Vector< std::string> &GetStatisticsTableHeaders();

      //! @brief function to calculate to statistics with regards to the given domain
      //! @param SSE_DOMAIN domain to be used for comparison against the given pool
      //! @return a table with a single row that stores the statistics regarding the pool
      storage::Table< double> CalculateStatistics( const DomainInterface &SSE_DOMAIN) const;

    private:

      //! @brief initialize from a vector of sses
      //! @param SSE_VECTOR vector of sses
      //! @param IGNORE_UNSTRUCTURED ignore unstructured sses
      //! @param IDEALIZE whether to idealize the SSEs (or simply move them to the origin)
      void Initialize
      (
        const util::SiPtrVector< const SSE> &SSE_VECTOR,
        const bool IGNORE_UNSTRUCTURED = true,
        const bool IDEALIZE = true
      );

    }; // class SSEPool

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_SSE_POOL_H_
