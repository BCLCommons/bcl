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
#include "assemble/bcl_assemble_sse_pool.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_stage_factory.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_range_set.h"
#include "pdb/bcl_pdb_line.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEPool::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPool())
    );

    //! @brief return command line flag for providing a pool file
    //! @return command line flag for providing a pool file
    util::ShPtr< command::FlagInterface> &SSEPool::GetFlagPoolRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "pool", "file with sse definitions to be used in constructing SSE pool",
          command::Parameter
          (
            "pool_filename", "filename for input SSE definitions to be used in constructing SSE pool", "pool.txt"
          )
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for setting minimum helix and strand lengths to be considered from the pool
    //! @return command line flag setting minimum helix and strand lengths to be considered from the pool
    util::ShPtr< command::FlagInterface> &SSEPool::GetFlagMinSSELengths()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "pool_min_sse_lengths", "minimal sse lengths for sses in the pool to be considered in folding"
        )
      );
      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_helix_length_param
      (
        new command::Parameter
        (
          "pool_min_helix_length", "\tminimal length of helices from the pool to be considered in folding", "9"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_strand_length_param
      (
        new command::Parameter
        (
          "pool_min_strand_length", "\tminimal length of strands from the pool to be considered in folding", "5"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters into flag
        flag->PushBack( s_helix_length_param);
        flag->PushBack( s_strand_length_param);
      }

      // end
      return s_flag;
    }

    //! @brief return command line flag for setting minimum helix and strand lengths to be considered from the pool
    //! @return command line flag setting minimum helix and strand lengths to be considered from the pool
    storage::Map< biol::SSType, size_t> SSEPool::GetCommandLineMinSSELengths()
    {
      // initialize sse min sizes for the pool
      storage::Map< biol::SSType, size_t> pool_sse_min_sizes;

      // set the values from the command line
      pool_sse_min_sizes[ biol::GetSSTypes().HELIX] =
          GetFlagMinSSELengths()->GetParameterList()( 0)->GetNumericalValue< size_t>();
      pool_sse_min_sizes[ biol::GetSSTypes().STRAND] =
          GetFlagMinSSELengths()->GetParameterList()( 1)->GetNumericalValue< size_t>();

      // end
      return pool_sse_min_sizes;
    }

    //! @brief return command line flag for specifying a pool prefix - to be used when different pools postfixes are
    //!        specified within a stages file
    //! @return command line flag setting the prefix of pool files that will be used
    util::ShPtr< command::FlagInterface> &SSEPool::GetFlagPoolPrefix()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "pool_prefix", "return command line flag for specifying a pool prefix - to be used when different pools "
          " postfixes are specified within a stages file. This flag does nothing unless the "
          + util::Format()( fold::StageFactory::e_PoolPostfix) + " tag is used in a stages file. In which case, the "
          " pool file will be "
          "searched for as \"<pool_prefix>.<pool_postfix_from_stages_file>\". The period (\".\") will be added "
          "automatically so do not include this in either the prefix or postfix."
        )
      );
      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_param
      (
        new command::Parameter
        (
          "pool_prefix", "\tstring which is the prefix for accessing the desired pool", "myproteinpoolprefix"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters into flag
        flag->PushBack( s_param);
      }

      // end
      return s_flag;

    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPool::SSEPool() :
      m_Data()
    {
    }

    //! @brief constructor from a vector of sses
    //! @param SSE_VECTOR vector of sses
    //! @param IGNORE_UNSTRUCTURED ignore unstructured sses
    //! @param IDEALIZE whether to idealize the SSEs (or simply move them to the origin)
    SSEPool::SSEPool
    (
      const util::SiPtrVector< const SSE> &SSE_VECTOR,
      const bool IGNORE_UNSTRUCTURED,
      const bool IDEALIZE
    ) :
      m_Data()
    {
      Initialize( SSE_VECTOR, IGNORE_UNSTRUCTURED, IDEALIZE);
    }

    //! @brief constructor from a list of sses
    //! @param SSE_LIST list of sses
    //! @param IGNORE_UNSTRUCTURED ignore unstructured sses
    //! @param IDEALIZE whether to idealize the SSEs (or simply move them to the origin)
    SSEPool::SSEPool
    (
      const util::SiPtrList< const SSE> &SSE_LIST,
      const bool IGNORE_UNSTRUCTURED,
      const bool IDEALIZE
    ) :
      m_Data()
    {
      Initialize( util::SiPtrVector< const SSE>( SSE_LIST.Begin(), SSE_LIST.End()), IGNORE_UNSTRUCTURED, IDEALIZE);
    }

    //! @brief virtual copy constructor
    SSEPool *SSEPool::Clone() const
    {
      return new SSEPool( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPool::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns SiPtrVector of all SSEs
    //! @return SiPtrVector of all SSEs
    util::SiPtrVector< const SSE> SSEPool::GetSSEs() const
    {
      util::SiPtrVector< const SSE> sses;

      //loop over all SSElements
      for
      (
        const_iterator set_itr( m_Data.Begin()), set_itr_end( m_Data.End());
        set_itr != set_itr_end;
        ++set_itr
      )
      {
        sses.PushBack( &( ( **set_itr)));
      }

      return sses;
    }

    //! @brief Get the total number of potentially structured AAs
    size_t SSEPool::GetNumberPotentiallyStructuredAAs() const
    {
      // set of structured ranges
      math::RangeSet< int> structured_ranges;
      //loop over all SSElements
      for
      (
        const_iterator set_itr( m_Data.Begin()), set_itr_end( m_Data.End());
        set_itr != set_itr_end;
        ++set_itr
      )
      {
        structured_ranges +=
          math::Range< int>( ( ( **set_itr).GetFirstAA())->GetSeqID(), ( ( **set_itr).GetLastAA())->GetSeqID());
      }
      size_t number_structured( 0);
      for
      (
        storage::Set< math::Range< int> >::const_iterator
          itr( structured_ranges.GetRanges().Begin()), itr_end( structured_ranges.GetRanges().End());
        itr != itr_end;
        ++itr
      )
      {
        // offset of 1 because the range borders are closed, so the last position would otherwise not be counted
        number_structured += itr->GetWidth() + 1;
      }
      return number_structured;
    }

    //! @brief returns SiPtrVector of all SSEs of given SS_TYPE
    //! @param SS_TYPE SSType of interest
    //! @return SiPtrVector of all SSEs of given SS_TYPE
    util::SiPtrVector< const SSE> SSEPool::GetSSEs( const biol::SSType &SS_TYPE) const
    {
      util::SiPtrVector< const SSE> sses;

      //loop over all SSElements
      for
      (
        const_iterator set_itr( m_Data.Begin()), set_itr_end( m_Data.End());
        set_itr != set_itr_end;
        ++set_itr
      )
      {
        if( ( *set_itr)->GetType() == SS_TYPE)
        {
          sses.PushBack( ( **set_itr));
        }
      }

      return sses;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns whether the pool has overlapping SSEs
    //! @return whether the pool has overlapping SSEs
    bool SSEPool::IsOverlapping() const
    {
      // iterate over all SSEs in the pool
      for( const_iterator itr_a( m_Data.Begin()), itr_end( m_Data.End()); itr_a != itr_end; ++itr_a)
      {
        // create a second iterator
        const_iterator itr_b( itr_a);
        ++itr_b;

        for( ; itr_b != itr_end; ++itr_b)
        {
          // if overlapping return true
          if( biol::DoOverlap( **itr_a, **itr_b))
          {
            return true;
          }
        }
      }

      // if no overlaps found return false
      return false;
    }

    //! @brief removes the identical SSEs with the provided PROTEIN_MODEL from pool
    //! @param PROTEIN_MODEL protein model
    //! @return SiPtrList of SSEs in the pool excluding ones that are identical to SSEs in the provided PROTEIN_MODEL
    util::SiPtrList< const SSE> SSEPool::GetNonIdenticalSSEs
    (
      const ProteinModel &PROTEIN_MODEL
    ) const
    {
      // get sses in the protein model and sort them in a set
      const util::SiPtrVector< const SSE> sses_in_model( PROTEIN_MODEL.GetSSEs());
      const storage::Set< util::SiPtr< const SSE>, SSELessThan> sses_in_model_sorted
      (
        sses_in_model.Begin(), sses_in_model.End()
      );

      // construct the list to be returned
      util::SiPtrList< const SSE> non_identical_sses;

      // take the difference of sses in pool and the protein model using the SSELessThan
      std::set_difference
      (
        m_Data.Begin(),
        m_Data.End(),
        sses_in_model_sorted.Begin(),
        sses_in_model_sorted.End(),
        std::inserter( non_identical_sses.InternalData(), non_identical_sses.InternalData().begin()),
        SSELessThan()
      );

      // return non_identical_sses
      return non_identical_sses;
    }

    //! @brief removes the overlapping SSEs with the provided PROTEIN_MODEL from pool
    //! @param PROTEIN_MODEL protein model
    //! @return SiPtrList of SSEs where there are no overlapping SSEs with the provided PROTEIN_MODEL
    util::SiPtrList< const SSE> SSEPool::GetNonOverlappingSSEs
    (
      const ProteinModel &PROTEIN_MODEL
    ) const
    {

      //TODO get Chain directly and do a find instead of iteration
      // initialize the list to be returned
      util::SiPtrList< const SSE> sse_list;

      // get the sses in the protein model
      util::SiPtrVector< const SSE> sses_from_model( PROTEIN_MODEL.GetSSEs());

      // iterate over sses in this pool
      for
      (
        const_iterator sse_itr_a( m_Data.Begin()),
          sse_itr_end_a( m_Data.End());
        sse_itr_a != sse_itr_end_a;
        ++sse_itr_a
      )
      {
        // initialize bool
        bool overlaps( false);

        // iterate over the sses in the protein model
        for
        (
          util::SiPtrVector< const SSE>::const_iterator sse_itr_b( sses_from_model.Begin()),
            sse_itr_end_b( sses_from_model.End());
          sse_itr_b != sse_itr_end_b;
          ++sse_itr_b
        )
        {
          // if these two sses overlap set the flag and break
          if( biol::DoOverlap( **sse_itr_a, **sse_itr_b))
          {
            overlaps = true; break;
          }
        }

        // if the sse did not overlap with any sse from the protein model
        if( !overlaps)
        {
          // insert this sse into the pool to be returned
          sse_list.PushBack( *sse_itr_a);
        }

      }

      // end
      return sse_list;
    }

    //! @brief return the SSEs of the given type in the pool which do not overlap with SSEs in the given model
    //! @param PROTEIN_MODEL protein model to compare to
    //! @param SSE_TYPE which SSE type to consider
    //! @return SSEs of the given type in the pool which do not overlap with SSEs in the given model
    util::SiPtrList< const SSE> SSEPool::GetNonOverlappingSSEs
    (
      const ProteinModel &PROTEIN_MODEL,
      const biol::SSType &SSE_TYPE
    ) const
    {
      //TODO get Chain directly and do a find instead of iteration
       // initialize the list to be returned
       util::SiPtrList< const SSE> sse_list;

       // get the sses in the protein model
       util::SiPtrVector< const SSE> sses_from_model( PROTEIN_MODEL.GetSSEs());

       // iterate over sses in this pool
       for
       (
         const_iterator sse_itr_a( m_Data.Begin()),
           sse_itr_end_a( m_Data.End());
         sse_itr_a != sse_itr_end_a;
         ++sse_itr_a
       )
       {
         // initialize bool
         bool overlaps( false);

         // iterate over the sses in the protein model
         for
         (
           util::SiPtrVector< const SSE>::const_iterator sse_itr_b( sses_from_model.Begin()),
             sse_itr_end_b( sses_from_model.End());
           sse_itr_b != sse_itr_end_b;
           ++sse_itr_b
         )
         {
           // if these two sses overlap set the flag and break
           if( biol::DoOverlap( **sse_itr_a, **sse_itr_b))
           {
             overlaps = true; break;
           }
         }

         // if the sse did not overlap with any sse from the protein model and is of the given type
         if( !overlaps && ( **sse_itr_a).GetType() == SSE_TYPE)
         {
           // insert this sse into the pool to be returned
           sse_list.PushBack( *sse_itr_a);
         }
       }

       return sse_list;
    }

    //! @brief returns overlapping SSEs in the pool with the SSEs from provided SSE_DOMAIN
    //! @param SSE_DOMAIN domain
    //! @param EXCLUDE_IDENTICAL boolean to decide whether SSEs from domain should be included in the return
    //! @param EXCLUDE_DIFFERENT_SSTYPE whether to exclude overlapping SSEs that differ in SSType
    //! @return SiPtrList of SSEs which overlap with SSEs with the provided SSE_DOMAIN
    util::SiPtrList< const SSE> SSEPool::GetOverlappingSSEs
    (
      const DomainInterface &SSE_DOMAIN,
      const bool EXCLUDE_IDENTICAL,
      const bool EXCLUDE_DIFFERENT_SSTYPE
    ) const
    {
      // initialize the list to be returned
      util::SiPtrList< const SSE> sse_list;

      // get the sses in the protein model
      util::SiPtrVector< const SSE> sses_from_model( SSE_DOMAIN.GetSSEs());

      // iterate over sses in this pool
      for
      (
        const_iterator sse_itr_a( m_Data.Begin()),
          sse_itr_end_a( m_Data.End());
        sse_itr_a != sse_itr_end_a;
        ++sse_itr_a
      )
      {
        // iterate over the sses in the protein model
        for
        (
          util::SiPtrVector< const SSE>::const_iterator sse_itr_b( sses_from_model.Begin()),
            sse_itr_end_b( sses_from_model.End());
          sse_itr_b != sse_itr_end_b;
          ++sse_itr_b
        )
        {
          // if these two sses overlap then
          if( biol::DoOverlap( **sse_itr_a, **sse_itr_b))
          {
            // if the EXCLUDE_IDENTICAL flag is set and they are same continue
            if( EXCLUDE_IDENTICAL && **sse_itr_a == **sse_itr_b)
            {
              continue;
            }

            // if exclude different SSType is given and the types actually differ
            if( EXCLUDE_DIFFERENT_SSTYPE && ( *sse_itr_a)->GetType() != ( *sse_itr_b)->GetType())
            {
              continue;
            }

            // insert this sse to sse_list
            sse_list.PushBack( *sse_itr_a);

            // and break
            break;
          }
        }
      }

      // end
      return sse_list;
    }

    //! @brief returns overlapping SSEs in the pool with the provided TARGET_SSE
    //! @param TARGET_SSE SSE of interest
    //! @param EXCLUDE_IDENTICAL boolean to decide whether SSEs from ProteinModel should be included in the return
    //! @param EXCLUDE_DIFFERENT_SSTYPE whether to exclude overlapping SSEs that differ in SSType
    //! @return SiPtrList of SSEs which overlap with TARGET_SSE
    util::SiPtrList< const SSE> SSEPool::GetOverlappingSSEs
    (
      const SSE &TARGET_SSE,
      const bool EXCLUDE_IDENTICAL,
      const bool EXCLUDE_DIFFERENT_SSTYPE
    ) const
    {
      // initialize the list to be returned
      util::SiPtrList< const SSE> sse_list;

      // iterate over sses in this pool
      for
      (
        const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if exclude different SSType is given and the types actually differ
        if( EXCLUDE_DIFFERENT_SSTYPE && ( *sse_itr)->GetType() != TARGET_SSE.GetType())
        {
          continue;
        }

        // if the EXCLUDE_IDENTICAL flag is set and they are same continue
        if( EXCLUDE_IDENTICAL && **sse_itr == TARGET_SSE)
        {
          continue;
        }

        // if these two sses overlap then
        if( biol::DoOverlap( **sse_itr, TARGET_SSE))
        {
          // insert this sse to sse_list
          sse_list.PushBack( *sse_itr);
        }
      }

      // end
      return sse_list;
    }

    //! @brief converts the data in the set from SSELessThan to SSELessThanNoOverlap
    //! @return SSELessThanNoOverlap set
    storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> SSEPool::GetRandomNonOverlappingSet() const
    {
      // copy the initial SSEs in the set
      storage::Set< util::ShPtr< SSE>, SSELessThan> initial_set( m_Data);

      // initialize the set to be returned
      storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> non_overlapping_set;

      // loop until the initial set is empty
      while( !initial_set.IsEmpty())
      {
        // get a random iterator on the set
        storage::Set< util::ShPtr< SSE>, SSELessThan>::iterator sse_itr
        (
          random::GetGlobalRandom().Iterator( initial_set.Begin(), initial_set.End(), initial_set.GetSize())
        );

        // insert the sse into the non overlapping set (if possible)
        non_overlapping_set.Insert( **sse_itr);

        // remove the sse from the initial set
        initial_set.RemoveElement( sse_itr);
      }

      // end
      return non_overlapping_set;
    }

    //! @brief calculates the number of helices and strands in the pool
    //!        this will equal to actual counts if it's non overlapping pool
    //!        if overlapping, it will be average numbers that lead to a complete model with non-overlapping SSEs
    //! @return average number of helices and strands expected in a complete model from this pool
    storage::Pair< double, double> SSEPool::CalculateAverageHelixStrandCounts() const
    {
      // if the SSEs are not overlapping
      if( !IsOverlapping())
      {
        // then return the ratio directly
        return storage::Pair< double, double>
        (
          GetSSEs( biol::GetSSTypes().HELIX).GetSize(),
          GetSSEs( biol::GetSSTypes().STRAND).GetSize()
        );
      }

      // initialize counts
      storage::Pair< double, double> helix_strand_counts( 0.0, 0.0);

      // initialize the number of repeats
      const size_t nr_repeats( 10);

      // if overlapping then we have to get few random non-overlapping sets to calculate the expected average number
      for( size_t i = 0; i < nr_repeats; ++i)
      {
        // new overlapping set
        const storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> non_overlapping_set
        (
          GetRandomNonOverlappingSet()
        );

        // create a new non overlapping pool
        const SSEPool this_pool
        (
          util::SiPtrVector< const SSE>( non_overlapping_set.Begin(), non_overlapping_set.End())
        );

        // sum up the numbers
        helix_strand_counts.First() += ( this_pool.GetSSEs( biol::GetSSTypes().HELIX).GetSize());
        helix_strand_counts.Second() += ( this_pool.GetSSEs( biol::GetSSTypes().STRAND).GetSize());
      }

      // normalize the number
      helix_strand_counts.First() /= nr_repeats;
      helix_strand_counts.Second() /= nr_repeats;

      // end
      return helix_strand_counts;
    }

    //! @brief calculate an estimate of the helix to strand ratio and return it
    //! @return helix to strand ratio ( 0:1 if only strands, 1:0 if only helices)
    storage::Pair< double, double> SSEPool::CalculateHelixToStrandRatio() const
    {
      // get the number of strands and number of helices
      const double nr_helices( GetSSEs( biol::GetSSTypes().HELIX).GetSize());
      const double nr_strands( GetSSEs( biol::GetSSTypes().STRAND).GetSize());

      // make sure there are at least one helix or strand in the pool
      BCL_Assert
      (
        nr_helices + nr_strands > 0, "The pool needs to have at least one helix or strand to calculate ratio!"
      );

      // if only helices
      if( nr_helices > 0 && nr_strands == 0)
      {
        return storage::Pair< double, double>( 1.0, 0.0);
      }
      // if only strands
      if( nr_helices == 0 && nr_strands > 0)
      {
        return storage::Pair< double, double>( 0.0, 1.0);
      }

      // if the loop is not overlapping
      if( !IsOverlapping())
      {
        // then return the ratio directly
        return storage::Pair< double, double>( nr_helices / GetSize(), nr_strands / GetSize());
      }

      // get the average counts for helices and strands
      storage::Pair< double, double> helix_strand_counts( CalculateAverageHelixStrandCounts());

      // normalize by the total counts
      const double total_nr_sses( helix_strand_counts.First() + helix_strand_counts.Second());
      helix_strand_counts.First() /= total_nr_sses;
      helix_strand_counts.Second() /= total_nr_sses;

      // end
      return helix_strand_counts;
    }

    //! @brief number of sses above given length
    //! @param MIN_SIZE_HELIX
    //! @return nr strands size > MIN_SIZE_STRAND + nr helices size > MIN_HELIX_SIZE
    double SSEPool::CountLongNonOverlappingSSEs( const size_t MIN_SIZE_STRAND, const size_t MIN_SIZE_HELIX) const
    {
      double count( 0);

      const size_t nr_repeats( IsOverlapping() ? 10 : 1);

      // if overlapping then we have to get few random non-overlapping sets to calculate the expected average number
      for( size_t i = 0; i < nr_repeats; ++i)
      {
        // new overlapping set
        const storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> non_overlapping_set
        (
          GetRandomNonOverlappingSet()
        );

        // iterate over sses
        for
        (
            storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( non_overlapping_set.Begin()), sse_itr_end( non_overlapping_set.End());
          sse_itr != sse_itr_end;
          ++sse_itr
        )
        {
          if
          (
               ( ( *sse_itr)->GetType() == biol::GetSSTypes().STRAND && ( *sse_itr)->GetSize() > MIN_SIZE_STRAND)
            || ( ( *sse_itr)->GetType() == biol::GetSSTypes().HELIX && ( *sse_itr)->GetSize() > MIN_SIZE_HELIX)
          )
          {
            count += 1.0;
          }
        }
      }

      // end
      return count / nr_repeats;
    }

    //! @brief prune the pool from sses that are not in the given map or do not have the minimum length
    //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
    void SSEPool::Prune( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP)
    {
      // new set of sses
      storage::Set< util::ShPtr< SSE>, SSELessThan> new_data;

      // iterate over sses
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThan>::const_iterator sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        storage::Map< biol::SSType, size_t>::const_iterator min_itr( MIN_SSE_SIZE_MAP.Find( ( *sse_itr)->GetType()));

        // sstype in map?
        if( min_itr == MIN_SSE_SIZE_MAP.End())
        {
          continue;
        }

        // sse has minimal size
        if( ( *sse_itr)->GetSize() < min_itr->second)
        {
          continue;
        }

        // insert into new data
        new_data.Insert( *sse_itr);
      }

      // update data
      m_Data = new_data;
    }

    //! @brief chops the SSEs in of the type in the map into chunks of the minimal size given in the map
    //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
    void SSEPool::ChopSSEs( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP)
    {
      // make a copy of the pool
      SSEPool new_pool( *this);

      // iterate over SSEs in the pool
      for( const_iterator sse_itr( Begin()), sse_itr_end( End()); sse_itr != sse_itr_end; ++sse_itr)
      {
        // cannot chop if size is smaller than 3
        if( ( *sse_itr)->GetSize() < 3)
        {
          new_pool.Insert( *sse_itr);
          continue;
        }

        // check if that type is in the map
        const storage::Map< biol::SSType, size_t>::const_iterator ss_itr( MIN_SSE_SIZE_MAP.Find( ( *sse_itr)->GetType()));

        if( ss_itr == MIN_SSE_SIZE_MAP.End())
        {
          new_pool.Insert( *sse_itr);
          continue;
        }

        // store the set of chopped sses
        const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> chopped_sses
        (
          // chop sse into two parts that are at least m_MinSSELengths in length
          ( *sse_itr)->Chop( std::max( ss_itr->second, ( *sse_itr)->GetSize() / 2 - 1))
        );

        // chop up this sse and insert the resultant set of sses into SSEPool
        new_pool.InsertElements( chopped_sses.Begin(), chopped_sses.End());
      }

      // update data
      m_Data = new_pool.m_Data;
    }

    //! @brief join adjacent sses if the adjacent sse is shorter than the given sse min size and has the same type, starting from the outside of a sequence of adjacent sses
    //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
    //! @return true if at least one join was performed - one call might not perform all possible joins and it can be necessary to call that function multiple timess
    bool SSEPool::Join( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP)
    {
      // nothing to do for empty data
      if( m_Data.IsEmpty())
      {
        return false;
      }

      storage::Set< util::ShPtr< SSE>, SSELessThan> new_data;

      bool joined( false);

      // iterate over all sses
      for( const_iterator itr( m_Data.Begin()), itr_end( m_Data.End()); itr != itr_end;)
      {
        const SSE &current_sse( **itr);
        ++itr;

        // is the current sse long enough and should be joined
        const storage::Map< biol::SSType, size_t>::const_iterator find_itr( MIN_SSE_SIZE_MAP.Find( current_sse.GetType()));

        // if this was the last sse
        // or it is of type not to be considered for joining
        // or it is the sse before an sse with different type
        // or the following sse is not adjacent
        // insert to new data and move on
        if
        (
             itr == itr_end
          || find_itr == MIN_SSE_SIZE_MAP.End()
          || current_sse.GetType() != ( *itr)->GetType()
          || !current_sse.DoesPrecede( **itr)
        )
        {
          new_data.Insert( util::ShPtr< SSE>( current_sse.Clone()));
          continue;
        }

        // if either current or next sse is too short, merge them
        if( current_sse.GetSize() < find_itr->second || ( *itr)->GetSize() < find_itr->second)
        {
          joined = true;
          util::ShPtr< SSE> sp_sse( current_sse.Clone());
          sp_sse->AppendSequence( **itr, false);
          new_data.Insert( sp_sse);
          BCL_MessageDbg
          (
            "joining sses: " + current_sse.GetIdentification() + " " + ( *itr)->GetIdentification() +
            " into: " + sp_sse->GetIdentification()
          );
          ++itr;
        }
        else
        {
          new_data.Insert( util::ShPtr< SSE>( current_sse.Clone()));
        }
      }

      // update data
      m_Data = new_data;

      // return if something was joined
      return joined;
    }

    //! @brief separate adjacent sses of identical type, by 2 * nr residues, but only if the resulting sses still have the desired size
    //! @param MIN_SSE_SIZE_MAP map containing the desired ss types and their minimum length
    //! @param NR_RESIDUES number of residues
    //! @return true if at least one separation was performed
    bool SSEPool::Separate( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZE_MAP, const size_t NR_RESIDUES)
    {
      // nothing to do for empty data
      if( m_Data.IsEmpty() || NR_RESIDUES == 0)
      {
        return false;
      }

      storage::Set< util::ShPtr< SSE>, SSELessThan> new_data;
      const_iterator itr( m_Data.Begin()), itr_end( m_Data.End());
      util::ShPtr< SSE> sp_previous( *itr);
      ++itr;
      bool separated( false);

      // iterate over all sses
      for( ; itr != itr_end; ++itr)
      {
        util::ShPtr< SSE> sp_current( *itr);

        // are they of different type or do not precede each other
        if
        (
             !sp_previous->GetType()->IsStructured() || !sp_current->GetType()->IsStructured()
          || !sp_previous->DoesPrecede( *sp_current)
        )
        {
          // insert the previous, update the previous and continue
          new_data.Insert( sp_previous.HardCopy());
          sp_previous = sp_current;
          continue;
        }

        // sse type considered and the min size
        const storage::Map< biol::SSType, size_t>::const_iterator find_itr_prev( MIN_SSE_SIZE_MAP.Find( sp_previous->GetType()));
        const storage::Map< biol::SSType, size_t>::const_iterator find_itr_cur( MIN_SSE_SIZE_MAP.Find( sp_current->GetType()));

        if( find_itr_prev == MIN_SSE_SIZE_MAP.End() || find_itr_prev->second == 0)
        {
          // insert the previous and move on
          new_data.Insert( sp_previous.HardCopy());
          sp_previous = sp_current;
          continue;
        }

        biol::AASequence coil_seq;

        // sse type is considered
        // remove residues if remaining length permits
        if( sp_previous->GetSize() >= find_itr_prev->second + NR_RESIDUES)
        {
          separated = true;
          coil_seq.AppendSequence( sp_previous->SubSequence( sp_previous->GetSize() - NR_RESIDUES, NR_RESIDUES));
          sp_previous = util::ShPtr< SSE>( new SSE( sp_previous->SubSequence( 0, sp_previous->GetSize() - NR_RESIDUES), sp_previous->GetType()));
        }

        if( sp_current->GetSize() >= find_itr_cur->second + NR_RESIDUES)
        {
          separated = true;
          coil_seq.AppendSequence( sp_current->SubSequence( 0, NR_RESIDUES));
          sp_current = util::ShPtr< SSE>( new SSE( sp_current->SubSequence( NR_RESIDUES, sp_current->GetSize() - NR_RESIDUES), sp_current->GetType()));
        }

        // insert previous, coil and update
        new_data.Insert( sp_previous.HardCopy());
        if( coil_seq.GetSize() > 0)
        {
          new_data.Insert( util::ShPtr< SSE>( new SSE( coil_seq, biol::GetSSTypes().COIL)));
        }
        sp_previous = sp_current;
      }

      // insert the last element
      new_data.Insert( sp_previous.HardCopy());
      m_Data = new_data;

      // end
      return separated;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read SSEPool from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPool::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Data, ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief read SSEPool from std::istream from a specially formatted pool file
    //! @param ISTREAM input stream
    //! @param PROTEIN_MODEL ProteinModel which SSEs belong to
    //! @param MIN_HELIX_LENGTH Minimal length for a HELIX to be added to the pool
    //! @param MIN_STRAND_LENGTH Minimal length for a STRAND to be added to the pool
    //! @return istream which was read from
    std::istream &SSEPool::ReadSSEPool
    (
      std::istream &ISTREAM,
      const ProteinModel &PROTEIN_MODEL,
      const size_t MIN_HELIX_LENGTH,
      const size_t MIN_STRAND_LENGTH
    )
    {
      m_Data.Reset();

      // read the identifier
      util::ObjectInterface::ReadIdentifier( ISTREAM);

      // initialize buffer
      std::string buffer;

      // while reading lines from STREAM into buffer and EOF is not hit
      while( !ISTREAM.eof() && std::getline( ISTREAM, buffer))
      {
        // if buffer is empty break
        if( buffer.empty())
        {
          break;
        }
        pdb::LineType line_type( util::TrimString( buffer.substr( 0, 6)));

        // if the line type read is END, then break
        if( line_type == pdb::GetLineTypes().END)
        {
          break;
        }
        // initialize pdb line from buffer
        pdb::Line current_line( buffer);

        // if helix
        if( line_type == pdb::GetLineTypes().HELIX)
        {
          // find the chain with matching chain id
          const util::SiPtr< const Chain> this_chain
          (
            PROTEIN_MODEL.GetChain( current_line.GetChar( pdb::GetEntryTypes().HELIXChainID_Initial))
          );

          // assert the chain is found
          BCL_Assert
          (
            this_chain.IsDefined(),
            "Unable to find the chain with id " +
              util::Format()( current_line.GetChar( pdb::GetEntryTypes().HELIXChainID_Initial))
          );

          const size_t seq_id_start
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().HELIXSequenceID_Initial)
          );
          const size_t seq_id_end
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().HELIXSequenceID_Terminal)
          );

          // if this sse is not long enough skip
          if( seq_id_end - seq_id_start + 1 < MIN_HELIX_LENGTH)
          {
            continue;
          }

          // create sse from subsequence
          util::ShPtr< SSE> this_sse
          (
            new SSE
            (
              this_chain->GetSequence()->SubSequence
              (
                seq_id_start - 1,
                seq_id_end - seq_id_start + 1
              ),
              biol::GetSSTypes().HELIX
            )
          );

          // check that the seqids and names of sses in the pool and the created one is correct
          BCL_Assert
          (
            size_t( this_sse->GetFirstAA()->GetSeqID()) == seq_id_start &&
            size_t( this_sse->GetLastAA()->GetSeqID()) == seq_id_end &&
            this_sse->GetFirstAA()->GetType()->GetParentTypeString() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().HELIXResidueName_Initial)
              )->GetParentTypeString() &&
            this_sse->GetLastAA()->GetType()->GetParentTypeString() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().HELIXResidueName_Terminal)
              )->GetParentTypeString(),
            "The seqids and residue names do not match for pool and protein_model " + current_line.GetString() +
            " for SSE " + this_sse->GetType().GetName() + " " + this_sse->GetChainID() + "|" +
              util::Format()( this_sse->GetFirstAA()->GetSeqID()) + "|" +
              util::Format()( this_sse->GetFirstAA()->GetType()->GetParentTypeString()) + "|" +
              util::Format()( this_sse->GetLastAA()->GetSeqID()) + " " +
              util::Format()( this_sse->GetLastAA()->GetType()->GetParentTypeString())
          );

          // initialize sse coordinates to origin and idealize
          this_sse->SetToIdealConformationAtOrigin();

          // insert this sse into the set
          if( !m_Data.Insert( this_sse).second)
          {
            BCL_MessageCrt( "unable to insert sse:" + current_line.GetString());
          }
        }
        // else if strand
        else if( line_type == pdb::GetLineTypes().SHEET)
        {
          // find the chain with matching chain id
          const util::SiPtr< const Chain> this_chain
          (
            PROTEIN_MODEL.GetChain( current_line.GetChar( pdb::GetEntryTypes().SHEETChainID_Initial))
          );

          // assert the chain is found
          BCL_Assert
          (
            this_chain.IsDefined(),
            "Unable to find the chain with id " +
              util::Format()( current_line.GetChar( pdb::GetEntryTypes().SHEETChainID_Initial))
          );

          const size_t seq_id_start
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().SHEETSequenceID_Initial)
          );
          const size_t seq_id_end
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().SHEETSequenceID_Terminal)
          );

          // if this sse is not long enough skip
          if( seq_id_end - seq_id_start + 1 < MIN_STRAND_LENGTH)
          {
            continue;
          }

          // create sse from subsequence
          util::ShPtr< SSE> this_sse
          (
            new SSE
            (
              this_chain->GetSequence()->SubSequence
              (
                seq_id_start - 1,
                seq_id_end - seq_id_start + 1
              ),
              biol::GetSSTypes().STRAND
            )
          );

          // check that the seqids and names of sses in the pool and the created one is correct
//          BCL_Assert
//          (
//            size_t( this_sse->GetFirstAA()->GetSeqID()) == seq_id_start &&
//            size_t( this_sse->GetLastAA()->GetSeqID()) == seq_id_end &&
//            this_sse->GetFirstAA()->GetType()->GetParentTypeString() ==
//              biol::GetAATypes().AATypeFromThreeLetterCode
//              (
//                current_line.GetString( pdb::GetEntryTypes().SHEETResidueName_Initial)
//              )->GetParentTypeString() &&
//            this_sse->GetLastAA()->GetType()->GetParentTypeString() ==
//              biol::GetAATypes().AATypeFromThreeLetterCode
//              (
//                current_line.GetString( pdb::GetEntryTypes().SHEETResidueName_Terminal)
//              )->GetParentTypeString(),
//              "The seqids and residue names do not match for pool and protein_model " + current_line.GetString() +
//              " for SSE " + this_sse->GetType().GetName() + " " + this_sse->GetChainID() + " " +
//                util::Format()( this_sse->GetFirstAA()->GetSeqID()) + " " +
//                util::Format()( this_sse->GetFirstAA()->GetType()->GetParentTypeString()) + " " +
//                util::Format()( this_sse->GetLastAA()->GetSeqID()) + " " +
//                util::Format()( this_sse->GetLastAA()->GetType()->GetParentTypeString())
//          );
          BCL_Assert
          (
            size_t( this_sse->GetFirstAA()->GetSeqID()) == seq_id_start,
              "The seqids start end " + current_line.GetString() +
              " for SSE " + this_sse->GetType().GetName() + " " + this_sse->GetChainID() + " " +
                util::Format()( this_sse->GetFirstAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetFirstAA()->GetType()->GetParentTypeString()) + " " +
                util::Format()( this_sse->GetLastAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetLastAA()->GetType()->GetParentTypeString())
          );

          BCL_Assert
          (

            size_t( this_sse->GetLastAA()->GetSeqID()) == seq_id_end
            ,
              "The seqids end " + current_line.GetString() +
              " for SSE " + this_sse->GetType().GetName() + " " + this_sse->GetChainID() + " " +
                util::Format()( this_sse->GetFirstAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetFirstAA()->GetType()->GetParentTypeString()) + " " +
                util::Format()( this_sse->GetLastAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetLastAA()->GetType()->GetParentTypeString())
          );
          BCL_Assert
          (

            this_sse->GetFirstAA()->GetType()->GetParentTypeString() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().SHEETResidueName_Initial)
              )->GetParentTypeString()
            ,
              "The aa type strings " + current_line.GetString() +
              " for SSE " + this_sse->GetType().GetName() + " " + this_sse->GetChainID() + " " +
                util::Format()( this_sse->GetFirstAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetFirstAA()->GetType()->GetParentTypeString()) + " " +
                util::Format()( this_sse->GetLastAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetLastAA()->GetType()->GetParentTypeString())
          );

          BCL_Assert
          (
            this_sse->GetLastAA()->GetType()->GetParentTypeString() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().SHEETResidueName_Terminal)
              )->GetParentTypeString(),
              "The aa type last " + current_line.GetString() + "|" + biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().SHEETResidueName_Terminal)
              )->GetParentTypeString() +
              "| for SSE " + this_sse->GetType().GetName() + " " + this_sse->GetChainID() + " " +
                util::Format()( this_sse->GetFirstAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetFirstAA()->GetType()->GetParentTypeString()) + " " +
                util::Format()( this_sse->GetLastAA()->GetSeqID()) + " " +
                util::Format()( this_sse->GetLastAA()->GetType()->GetParentTypeString())
          );
          // initialize sse coordinates to origin and idealize
          this_sse->SetToIdealConformationAtOrigin();

          // insert this sse into the set
          if( !m_Data.Insert( this_sse).second)
          {
            BCL_MessageCrt( "unable to insert sse:" + current_line.GetString());
          }
        }

        // if an unrecognized line type is supplied, break
        else
        {
          BCL_MessageCrt
          (
            "The following line type is not a helix or sheet " + line_type.GetName()
          );
          break;
        };
      } // while stream

      // end
      return ISTREAM;
    }

    //! @brief read SSEPool from std::istream from a specially formatted pool file
    //! @param ISTREAM input stream
    //! @param SEQUENCE aa sequence to which SSEs belong to
    //! @param MIN_HELIX_LENGTH Minimal length for a HELIX to be added to the pool
    //! @param MIN_STRAND_LENGTH Minimal length for a STRAND to be added to the pool
    //! @return istream which was read from
    std::istream &SSEPool::ReadSSEPool
    (
      std::istream &ISTREAM,
      const biol::AASequence &SEQUENCE,
      const size_t MIN_HELIX_LENGTH,
      const size_t MIN_STRAND_LENGTH
    )
    {
      // read the identifier
      util::ObjectInterface::ReadIdentifier( ISTREAM);

      // initialize buffer
      std::string buffer;

      // while reading lines from STREAM into buffer and EOF is not hit
      while( !ISTREAM.eof() && std::getline( ISTREAM, buffer))
      {
        // if buffer is empty break
        if( buffer.empty())
        {
          break;
        }
        pdb::LineType line_type( util::TrimString( buffer.substr( 0, 6)));

        // if the line type read is END, then break
        if( line_type == pdb::GetLineTypes().END)
        {
          break;
        }
        // initialize pdb line from buffer
        pdb::Line current_line( buffer);

        // if helix
        if( line_type == pdb::GetLineTypes().HELIX)
        {
          // skip if sequence chain id does not agree
          if( SEQUENCE.GetChainID() != current_line.GetChar( pdb::GetEntryTypes().HELIXChainID_Initial))
          {
            continue;
          }

          const size_t seq_id_start
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().HELIXSequenceID_Initial)
          );
          const size_t seq_id_end
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().HELIXSequenceID_Terminal)
          );

          // if this sse is not long enough skip
          if( seq_id_end - seq_id_start + 1 < MIN_HELIX_LENGTH)
          {
            continue;
          }

          // create sse from subsequence
          util::ShPtr< SSE> this_sse
          (
            new SSE
            (
              SEQUENCE.SubSequence
              (
                seq_id_start - 1,
                seq_id_end - seq_id_start + 1
              ),
              biol::GetSSTypes().HELIX
            )
          );

          // check that the seqids and names of sses in the pool and the created one is correct
          BCL_Assert
          (
            size_t( this_sse->GetFirstAA()->GetSeqID()) == seq_id_start &&
            size_t( this_sse->GetLastAA()->GetSeqID()) == seq_id_end &&
            this_sse->GetFirstAA()->GetType()->GetParentType() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().HELIXResidueName_Initial)
              )->GetParentType() &&
            this_sse->GetLastAA()->GetType()->GetParentType() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().HELIXResidueName_Terminal)
              )->GetParentType(),
            "The seqids and residue names do not match for pool and protein_model " + current_line.GetString() +
            " for SSE " + this_sse->GetIdentification()
          );

          // initialize sse coordinates to origin and idealize
          this_sse->SetToIdealConformationAtOrigin();

          // insert this sse into the set
          if( !m_Data.Insert( this_sse).second)
          {
            BCL_MessageCrt( "unable to insert sse:" + current_line.GetString());
          }
        }
        // else if strand
        else if( line_type == pdb::GetLineTypes().SHEET)
        {
          // find the chain with matching chain id
          if( SEQUENCE.GetChainID() != current_line.GetChar( pdb::GetEntryTypes().SHEETChainID_Initial))
          {
            continue;
          }

          const size_t seq_id_start
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().SHEETSequenceID_Initial)
          );
          const size_t seq_id_end
          (
            current_line.GetNumericalValue< size_t>( pdb::GetEntryTypes().SHEETSequenceID_Terminal)
          );

          // if this sse is not long enough skip
          if( seq_id_end - seq_id_start + 1 < MIN_STRAND_LENGTH)
          {
            continue;
          }

          // create sse from subsequence
          util::ShPtr< SSE> this_sse
          (
            new SSE
            (
              SEQUENCE.SubSequence
              (
                seq_id_start - 1,
                seq_id_end - seq_id_start + 1
              ),
              biol::GetSSTypes().STRAND
            )
          );

          // check that the seqids and names of sses in the pool and the created one is correct
          BCL_Assert
          (
            size_t( this_sse->GetFirstAA()->GetSeqID()) == seq_id_start &&
            size_t( this_sse->GetLastAA()->GetSeqID()) == seq_id_end &&
            this_sse->GetFirstAA()->GetType()->GetParentType() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().SHEETResidueName_Initial)
              )->GetParentType() &&
            this_sse->GetLastAA()->GetType()->GetParentType() ==
              biol::GetAATypes().AATypeFromThreeLetterCode
              (
                current_line.GetString( pdb::GetEntryTypes().SHEETResidueName_Terminal)
              )->GetParentType(),
              "The seqids and residue names do not match for pool and protein_model " + current_line.GetString() +
              " for SSE " + this_sse->GetIdentification()
          );

          // initialize sse coordinates to origin and idealize
          this_sse->SetToIdealConformationAtOrigin();

          // insert this sse into the set
          if( !m_Data.Insert( this_sse).second)
          {
            BCL_MessageCrt( "unable to insert sse:" + current_line.GetString());
          }
        }

        // if an unrecognized line type is supplied, break
        else
        {
          BCL_MessageCrt
          (
            "The following line type is not a helix or sheet " + line_type.GetName()
          );
          break;
        };
      } // while stream

      // end
      return ISTREAM;
    }

    //! @brief ReadSSEPoolInformation reads the pool file and gives the information it contains
    //!        This function does not create an SSE pool.
    //! @param POOL_FILENAME is the name and path of the file which contains the SSE pool
    //! @return returns a storage::Vector of the pdb::Lines which are contained in the SSE pool
    storage::Vector< pdb::Line> SSEPool::ReadSSEPoolInformation
    (
      const std::string &POOL_FILENAME
    )
    {
      // create io::IFStream "read" to be used to read in "POOL_FILENAME"
      io::IFStream read;

      // open "read" and bind it to "POOL_FILENAME"
      io::File::MustOpenIFStream( read, POOL_FILENAME.c_str());

      // read the identifier
      SSEPool().ReadIdentifier( read);

      // create std::string "buffer" which will be used to temporarily hold each line as it is read in
      std::string buffer;

      // create storage::Vector< pdb::Line> "pdb_lines" which will be used to hold all the pdb::Lines in the pool file
      storage::Vector< pdb::Line> pdb_lines;

      // get lines from "read" and put into "buffer" while EOF is not hit and a line can be gotten
      while( !read.eof() && std::getline( read, buffer))
      {
        // if buffer is empty break
        if( buffer.empty())
        {
          break;
        }

        // insert the current pool line contained in "buffer" into "pdb_lines"
        pdb_lines.PushBack( pdb::Line( buffer));
      }

      // close and clear the stream
      io::File::CloseClearFStream( read);

      // return the vector of pdb lines
      return pdb_lines;
    }

    //! @brief GetChainsRepresented obtains a storage::Set of all the chain ids in a pool file
    //!        This function does not create an SSE pool.
    //! @param POOL_FILENAME is the name and path of the file which contains the SSE pool
    //! @return returns a storage::Set< char> which has all the chain ids contained in the SSE pool
    storage::Set< char> SSEPool::GetChainsRepresented
    (
      const std::string &POOL_FILENAME
    )
    {
      // create vector of pdb::Lines "pdb_lines" and initialize with the lines in "POOL_FILENAME"
      const storage::Vector< pdb::Line> pdb_lines( ReadSSEPoolInformation( POOL_FILENAME));

      // create Set "chain_ids" to hold the chain ids that are in "pdb_lines"
      storage::Set< char> chain_ids;

      // iterate through "pdb_lines" to get the chain ids
      for
      (
        storage::Vector< pdb::Line>::const_iterator line_itr( pdb_lines.Begin()), line_itr_end( pdb_lines.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        // true if pdb::Line currently denoted by "line_itr" is of type HELIX
        if( line_itr->GetType() == pdb::GetLineTypes().HELIX)
        {
          // create const char "helix_id_initial" and "helix_id_terminal" initialize with chain ids from "line_itr"
          const char helix_id_initial( line_itr->GetChar( pdb::GetEntryTypes().HELIXChainID_Initial)),
                     helix_id_terminal( line_itr->GetChar( pdb::GetEntryTypes().HELIXChainID_Terminal));

          // make sure that the helix is contained within a single chain
          BCL_Assert
          (
            helix_id_initial == helix_id_terminal,
            "Initial helix id " + util::Format()( helix_id_initial) + " does not match terminal helix id " +
            util::Format()( helix_id_terminal)
          );

          // insert "strand_id_initial" into "chain_ids"
          chain_ids.Insert( helix_id_initial);
        }
        // true if pdb::Line currently denoted by "line_itr" is of type STRAND
        else if( line_itr->GetType() == pdb::GetLineTypes().SHEET)
        {
          // create const char "strand_id_initial" and "strand_id_terminal" initialize with chain ids from "line_itr"
          const char strand_id_initial( line_itr->GetChar( pdb::GetEntryTypes().SHEETChainID_Initial)),
                     strand_id_terminal( line_itr->GetChar( pdb::GetEntryTypes().SHEETChainID_Terminal));

          // make sure that the sheet is contained within a single chain
          BCL_Assert
          (
            strand_id_initial == strand_id_terminal,
            "Initial strand id " + util::Format()( strand_id_initial) + " does not match terminal stand id " +
            util::Format()( strand_id_terminal)
          );

          // insert "strand_id_initial" into "chain_ids"
          chain_ids.Insert( strand_id_initial);
        }
      }

      // return "chain_ids" which has the chain ids in "pdb_lines"
      return chain_ids;
    }

    //! @brief write to SSEPool std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &SSEPool::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

    //! @brief write SSEPool to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &SSEPool::WriteSSEPool( std::ostream &OSTREAM) const
    {
      // output identifier first
      util::ObjectInterface::WriteIdentifier( OSTREAM, 0);

      // serial of sses when written
      size_t sse_serial( 1);

      // iterate over the pool
      for
      (
        const_iterator sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // Write to pdb line
        const util::ShPtr< pdb::Line> this_line( pdb::Factory::WriteSSEDefinitionToLine( **sse_itr, sse_serial));

        if( this_line.IsDefined())
        {
          // write to stream
          OSTREAM << this_line->GetString() << '\n';

          // increment sse serial
          ++sse_serial;
        }
        else
        {
          BCL_MessageDbg( "did not write SSE to pool: " + ( *sse_itr)->GetIdentification());
        }
      }

      // print out the end
      OSTREAM << pdb::Line( pdb::GetLineTypes().END).GetString() << '\n';

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief finds SSEs in pool that agree best with SSEs in PROTEIN_MODEL
    //! @param PROTEIN_MODEL protein model whose SSEs are to be checked against the pool
    //! @return ShPtrVector to SSEs in pool that agree best with SSEs in PROTEIN_MODEL
    util::SiPtrVector< const SSE> SSEPool::FindBestMatchFromPool( const ProteinModel &PROTEIN_MODEL) const
    {
      // initialize the SiPtrVector to be returned
      util::SiPtrVector< const SSE> best_matches_from_pool;

      // get the sses in the protein model
      util::SiPtrVector< const SSE> sses_from_model( PROTEIN_MODEL.GetSSEs());

      // iterate over the sses in the protein model
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr_protein_model( sses_from_model.Begin()),
          sse_itr_end_protein_model( sses_from_model.End());
        sse_itr_protein_model != sse_itr_end_protein_model;
        ++sse_itr_protein_model
      )
      {
        // store the SiPtr to matching SSE
        util::SiPtr< const SSE> this_match( FindBestMatchFromPool( **sse_itr_protein_model).First());

        // if match is defined
        if( this_match.IsDefined())
        {
          // pushback SiPtr to identified sses from pool into best_matches_from_pool
          best_matches_from_pool.PushBack( this_match);
        }
        else
        {
          BCL_MessageStd
          (
            "no appropriate sse in pool for real sse " +
            util::Format()( ( *sse_itr_protein_model)->GetFirstAA()->GetSeqID()) + " " +
            util::Format()( ( *sse_itr_protein_model)->GetLastAA()->GetSeqID())
          );
        }
      }

      // return SiPtrVector of best matches from pool
      return best_matches_from_pool;
    }

    //! @brief finds SSEs in pool that agree best with SSE that is passed to function
    //! @param TARGET_SSE SSE to be checked against the pool
    //! @param MATCH_THRESHOLD the multiplier for setting the minimum threshold to be considered a match. This value is
    //!                        multiplied by the length of the TARGET_SSE to get the minimum threshold value
    //! @return ShPtr to SSE in pool that agrees best with SSE that is passed to function
    storage::Pair< util::SiPtr< const SSE>, double>
    SSEPool::FindBestMatchFromPool( const SSE &TARGET_SSE, const double MATCH_THRESHOLD) const
    {
      // initialize SiPtr to be returned
      storage::Pair< util::SiPtr< const SSE>, double> best_match_from_pool;
      best_match_from_pool.Second() = MATCH_THRESHOLD * double( TARGET_SSE.GetSize());

      // get list of SSEs in pool that overlap with TARGET_SSE
      util::SiPtrList< const SSE> overlapping_sse_list( GetOverlappingSSEs( TARGET_SSE, false));

      if( overlapping_sse_list.IsEmpty())
      {
        BCL_MessageDbg( "GetOverlappingSSEs function didn't find any match ");
      }

      // iterate over sses in overlapping_sse_list
      for
      (
        util::SiPtrList< const SSE>::const_iterator sse_itr( overlapping_sse_list.Begin()),
          sse_itr_end( overlapping_sse_list.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // skip all the sse pairs that are not of the same type
        if( TARGET_SSE.GetType() != ( *sse_itr)->GetType())
        {
          continue;
        }

        // calculate "temp overlap measure" between TARGET_SSE and sse from list
        const double temp_overlap_measure( CalculateOverlapMeasure( TARGET_SSE, **sse_itr));

        if( temp_overlap_measure < best_match_from_pool.Second())
        {
          best_match_from_pool.First() = *sse_itr;
          best_match_from_pool.Second() = temp_overlap_measure;
        }
      }

      BCL_MessageDbg
      (
        "best_overlap_measure " + util::Format()( best_match_from_pool.Second()) +
         " , " + util::Format()( TARGET_SSE.GetType().GetName())
      );

      return best_match_from_pool;
    }

    //! @brief calculates overlap measure between TARGET_SSE and SSE_FROM_POOL
    //! @param TARGET_SSE SSE to be checked against SSE from pool
    //! @param SSE_FROM_POOL SSE to be checked against TARGET_SSE
    //! @return int overlap measure
    int SSEPool::CalculateOverlapMeasure( const SSE &TARGET_SSE, const SSE &SSE_FROM_POOL) const
    {
      // total number of residues that the two SSEs differ in
      int residue_difference
      (
        math::Absolute( TARGET_SSE.GetFirstAA()->GetSeqID() - SSE_FROM_POOL.GetFirstAA()->GetSeqID()) +
        math::Absolute( TARGET_SSE.GetLastAA()->GetSeqID() - SSE_FROM_POOL.GetLastAA()->GetSeqID())
      );

      // difference in length between the two SSEs
      int length_difference( math::Absolute( int( TARGET_SSE.GetSize()) - int( SSE_FROM_POOL.GetSize())));

      // calculate overlap measure as sum residue_difference and length_difference
      int overlap_measure( residue_difference + length_difference);

      return overlap_measure;
    }

    //! @brief function to return statistics table header names
    //! @return vector of table header names
    const storage::Vector< std::string> &SSEPool::GetStatisticsTableHeaders()
    {
      // initialize static variable
      static storage::Vector< std::string> header_names;

      // if not initialized yet
      if( header_names.IsEmpty())
      {
        header_names.PushBack( "nr_sse");
        header_names.PushBack( "nr_helix");
        header_names.PushBack( "nr_strand");
        header_names.PushBack( "nr_sse_pool");
        header_names.PushBack( "nr_helix_pool");
        header_names.PushBack( "nr_strand_pool");
        header_names.PushBack( "nr_identified_sses");
        header_names.PushBack( "nr_identified_helix");
        header_names.PushBack( "nr_identified_strand");
        header_names.PushBack( "nr_missing_sses");
        header_names.PushBack( "nr_missing_helix");
        header_names.PushBack( "nr_missing_strand");
        header_names.PushBack( "nr_correct_sses_pool");
        header_names.PushBack( "nr_correct_helix_pool");
        header_names.PushBack( "nr_correct_strand_pool");
        header_names.PushBack( "nr_false_sses_pool");
        header_names.PushBack( "nr_false_helix_pool");
        header_names.PushBack( "nr_false_strand_pool");
        header_names.PushBack( "nr_overlaps");
        header_names.PushBack( "nr_helix_overlaps");
        header_names.PushBack( "nr_strand_overlaps");
        header_names.PushBack( "best_avg_overlap");
        header_names.PushBack( "best_avg_helix_overlap");
        header_names.PushBack( "best_avg_strand_overlap");
        header_names.PushBack( "best_avg_shift");
        header_names.PushBack( "best_avg_helix_shift");
        header_names.PushBack( "best_avg_strand_shift");
        header_names.PushBack( "best_avg_length_dev");
        header_names.PushBack( "best_avg_helix_length_dev");
        header_names.PushBack( "best_avg_strand_length_dev");
        header_names.PushBack( "all_avg_overlap");
        header_names.PushBack( "all_avg_helix_overlap");
        header_names.PushBack( "all_avg_strand_overlap");
        header_names.PushBack( "all_avg_shift");
        header_names.PushBack( "all_avg_helix_shift");
        header_names.PushBack( "all_avg_strand_shift");
        header_names.PushBack( "all_avg_length_dev");
        header_names.PushBack( "all_avg_helix_length_dev");
        header_names.PushBack( "all_avg_strand_length_dev");
        header_names.PushBack( "best_sum_overlap");
        header_names.PushBack( "best_sum_helix_overlap");
        header_names.PushBack( "best_sum_strand_overlap");
        header_names.PushBack( "best_sum_shift");
        header_names.PushBack( "best_sum_helix_shift");
        header_names.PushBack( "best_sum_strand_shift");
        header_names.PushBack( "best_sum_length_dev");
        header_names.PushBack( "best_sum_helix_length_dev");
        header_names.PushBack( "best_sum_strand_length_dev");
        header_names.PushBack( "all_sum_overlap");
        header_names.PushBack( "all_sum_helix_overlap");
        header_names.PushBack( "all_sum_strand_overlap");
        header_names.PushBack( "all_sum_shift");
        header_names.PushBack( "all_sum_helix_shift");
        header_names.PushBack( "all_sum_strand_shift");
        header_names.PushBack( "all_sum_length_dev");
        header_names.PushBack( "all_sum_helix_length_dev");
        header_names.PushBack( "all_sum_strand_length_dev");
      }

      // end
      return header_names;
    }

    //! @brief function to calculate to statistics with regards to the given domain
    //! @param SSE_DOMAIN domain to be used for comparison against the given pool
    //! @return a table with a single row that stores the statistics regarding the pool
    storage::Table< double> SSEPool::CalculateStatistics( const DomainInterface &SSE_DOMAIN) const
    {
      // construct table header
      const util::ShPtr< storage::TableHeader> sp_table_header
      (
        new storage::TableHeader( GetStatisticsTableHeaders())
      );

      // create table
      storage::Table< double> table( sp_table_header);

      // get the sses in the protein model
      util::SiPtrVector< const SSE> sses_from_model( SSE_DOMAIN.GetSSEs());

      // get SSEs from protein and construct another pool
      SSEPool native_pool( sses_from_model);

      // construct initial values
      const double nr_sse( native_pool.GetSize());
      const double nr_helix( native_pool.GetSSEs( biol::GetSSTypes().HELIX).GetSize());
      const double nr_strand( native_pool.GetSSEs( biol::GetSSTypes().STRAND).GetSize());
      const double nr_sse_pool( GetSize());
      const double nr_helix_pool( GetSSEs( biol::GetSSTypes().HELIX).GetSize());
      const double nr_strand_pool( GetSSEs( biol::GetSSTypes().STRAND).GetSize());

      // construct row_data with initial numbers
      storage::Vector< double> row_data
      (
        storage::Vector< double>::Create( nr_sse, nr_helix, nr_strand, nr_sse_pool, nr_helix_pool, nr_strand_pool)
      );

      // allocate memory
      row_data.AllocateMemory( sp_table_header->GetSize());

      // initialize sum for overlap
      double best_overlap_sum( 0.0);
      double best_overlap_helix_sum( 0.0);
      double best_overlap_strand_sum( 0.0);
      double best_shift_sum( 0.0);
      double best_shift_helix_sum( 0.0);
      double best_shift_strand_sum( 0.0);
      double best_length_dev_sum( 0.0);
      double best_length_dev_helix_sum( 0.0);
      double best_length_dev_strand_sum( 0.0);
      double all_overlap_sum( 0.0);
      double all_overlap_helix_sum( 0.0);
      double all_overlap_strand_sum( 0.0);
      double all_shift_sum( 0.0);
      double all_shift_helix_sum( 0.0);
      double all_shift_strand_sum( 0.0);
      double all_length_dev_sum( 0.0);
      double all_length_dev_helix_sum( 0.0);
      double all_length_dev_strand_sum( 0.0);

      //  initialize counts for overlaps
      double nr_overlaps( 0);
      double nr_helix_overlaps( 0);
      double nr_strand_overlaps( 0);

      // initialize counter for identified sses, helices and strands
      double nr_identified_sses( 0);
      double nr_identified_helix( 0);
      double nr_identified_strand( 0);

      // iterate over the sses in the protein model
      for
      (
        util::SiPtrVector< const SSE>::const_iterator
          sse_itr( sses_from_model.Begin()), sse_itr_end( sses_from_model.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // create reference on the SSE
        const SSE &this_sse( **sse_itr);

        BCL_MessageDbg( "Looking at SSE " + this_sse.GetIdentification());

        // initialize the values for the best overlap, shift and length dev for this SSE
        util::SiPtr< const SSE> best_match;
        double best_overlap( std::numeric_limits< double>::max());
        double best_shift( std::numeric_limits< double>::max());
        double best_length_dev( std::numeric_limits< double>::max());

        // get list of SSEs in pool that overlap with this SSE and exclude ones that do not match the type
        util::SiPtrList< const SSE> overlapping_sses( GetOverlappingSSEs( this_sse, false, true));

        // if no match
        if( overlapping_sses.IsEmpty())
        {
          continue;
          BCL_MessageDbg( "\tno overlapping SSEs");
        }

        // increment counters
        ++nr_identified_sses;
        // if helix
        if( this_sse.GetType() == biol::GetSSTypes().HELIX)
        {
          ++nr_identified_helix;
        }
        // else if strand
        else
        {
          ++nr_identified_strand;
        }

        // iterate over overlapping SSEs from pool
        for
        (
          util::SiPtrList< const SSE>::const_iterator
            pool_itr( overlapping_sses.Begin()), pool_itr_end( overlapping_sses.End());
          pool_itr != pool_itr_end; ++pool_itr
        )
        {
          // calculate the the shift, length deviation and the overlap measure
          const double this_shift
          (
            math::Absolute( this_sse.GetFirstAA()->GetSeqID() - ( *pool_itr)->GetFirstAA()->GetSeqID()) +
            math::Absolute( this_sse.GetLastAA()->GetSeqID() - ( *pool_itr)->GetLastAA()->GetSeqID())
          );
          const double this_length_dev( math::Absolute( int( this_sse.GetSize()) - int( ( *pool_itr)->GetSize())));
          const double this_overlap( this_shift + this_length_dev);

          BCL_MessageDbg
          (
            "\tvs " + ( *pool_itr)->GetIdentification() +
            "\toverlap:\t" + util::Format()( this_overlap) +
            "\tshift:\t" + util::Format()( this_shift) +
            "\tlength_dev:\t" + util::Format()( this_length_dev)
          );

          // if this overlap is better than the best one so far
          if( this_overlap < best_overlap)
          {
            best_match = *pool_itr;
            best_overlap = this_overlap;
            best_shift = this_shift;
            best_length_dev = this_length_dev;
          }

          // update the stats for all overlaps
          ++nr_overlaps;
          all_overlap_sum += this_overlap;
          all_shift_sum += this_shift;
          all_length_dev_sum += this_length_dev;

          // if helix
          if( ( *pool_itr)->GetType() == biol::GetSSTypes().HELIX)
          {
            all_overlap_helix_sum += this_overlap;
            all_length_dev_helix_sum += this_length_dev;
            all_shift_helix_sum += this_shift;
            ++nr_helix_overlaps;
          }
          else
          {
            all_overlap_strand_sum += this_overlap;
            all_length_dev_strand_sum += this_length_dev;
            all_shift_strand_sum += this_shift;
            ++nr_strand_overlaps;
          }
        }

        // assert a match was found
        BCL_Assert( best_match.IsDefined(), "A match should have been found");

        BCL_MessageDbg
        (
          "\tbest overlap " + best_match->GetIdentification() +
          "\toverlap:\t" + util::Format()( best_overlap) +
          "\tshift:\t" + util::Format()( best_shift) +
          "\tlength_dev:\t" + util::Format()( best_length_dev)
        );

        // update the sums
        best_overlap_sum += best_overlap;
        best_shift_sum += best_shift;
        best_length_dev_sum += best_length_dev;

        // if helix
        if( ( *sse_itr)->GetType() == biol::GetSSTypes().HELIX)
        {
          best_overlap_helix_sum += best_overlap;
          best_length_dev_helix_sum += best_length_dev;
          best_shift_helix_sum += best_shift;
        }
        else
        {
          best_overlap_strand_sum += best_overlap;
          best_length_dev_strand_sum += best_length_dev;
          best_shift_strand_sum += best_shift;
        }
      }

      // calculate the identified sse counts and add values to row data
      row_data.PushBack( nr_identified_sses);
      row_data.PushBack( nr_identified_helix);
      row_data.PushBack( nr_identified_strand);
      row_data.PushBack( nr_sse - nr_identified_sses);
      row_data.PushBack( nr_helix - nr_identified_helix);
      row_data.PushBack( nr_strand - nr_identified_strand);

      // get the overlapping SSEs from the pool
      const SSEPool overlapping_pool_sses( GetOverlappingSSEs( SSE_DOMAIN, false, true));

      // update counts for the correct and false sses from the pool
      const double nr_correct_sses_pool( overlapping_pool_sses.GetSize());
      const double nr_correct_helix_pool( overlapping_pool_sses.GetSSEs( biol::GetSSTypes().HELIX).GetSize());
      const double nr_correct_strand_pool( overlapping_pool_sses.GetSSEs( biol::GetSSTypes().STRAND).GetSize());
      row_data.PushBack( nr_correct_sses_pool);
      row_data.PushBack( nr_correct_helix_pool);
      row_data.PushBack( nr_correct_strand_pool);
      row_data.PushBack( nr_sse_pool - nr_correct_sses_pool);
      row_data.PushBack( nr_helix_pool - nr_correct_helix_pool);
      row_data.PushBack( nr_strand_pool - nr_correct_strand_pool);

      // counts for overlaps
      row_data.PushBack( nr_overlaps);
      row_data.PushBack( nr_helix_overlaps);
      row_data.PushBack( nr_strand_overlaps);

      // calculate and insert the average overlap measures for best overlaps
      row_data.PushBack( nr_identified_sses    == 0 ? 0.0 : best_overlap_sum / nr_identified_sses);
      row_data.PushBack( nr_identified_helix   == 0 ? 0.0 : best_overlap_helix_sum / nr_identified_helix);
      row_data.PushBack( nr_identified_strand  == 0 ? 0.0 : best_overlap_strand_sum / nr_identified_strand);
      row_data.PushBack( nr_identified_sses    == 0 ? 0.0 : best_shift_sum / nr_identified_sses);
      row_data.PushBack( nr_identified_helix   == 0 ? 0.0 : best_shift_helix_sum / nr_identified_helix);
      row_data.PushBack( nr_identified_strand  == 0 ? 0.0 : best_shift_strand_sum / nr_identified_strand);
      row_data.PushBack( nr_identified_sses    == 0 ? 0.0 : best_length_dev_sum / nr_identified_sses);
      row_data.PushBack( nr_identified_helix   == 0 ? 0.0 : best_length_dev_helix_sum / nr_identified_helix);
      row_data.PushBack( nr_identified_strand  == 0 ? 0.0 : best_length_dev_strand_sum / nr_identified_strand);

      // insert average overlap measures for all overlaps
      row_data.PushBack( nr_overlaps        == 0 ? 0.0 : all_overlap_sum / nr_overlaps);
      row_data.PushBack( nr_helix_overlaps  == 0 ? 0.0 : all_overlap_helix_sum / nr_helix_overlaps);
      row_data.PushBack( nr_strand_overlaps == 0 ? 0.0 : all_overlap_strand_sum / nr_strand_overlaps);
      row_data.PushBack( nr_overlaps        == 0 ? 0.0 : all_shift_sum / nr_overlaps);
      row_data.PushBack( nr_helix_overlaps  == 0 ? 0.0 : all_shift_helix_sum / nr_helix_overlaps);
      row_data.PushBack( nr_strand_overlaps == 0 ? 0.0 : all_shift_strand_sum / nr_strand_overlaps);
      row_data.PushBack( nr_overlaps        == 0 ? 0.0 : all_length_dev_sum / nr_overlaps);
      row_data.PushBack( nr_helix_overlaps  == 0 ? 0.0 : all_length_dev_helix_sum / nr_helix_overlaps);
      row_data.PushBack( nr_strand_overlaps == 0 ? 0.0 : all_length_dev_strand_sum / nr_strand_overlaps);

      // insert the sums
      row_data.PushBack( best_overlap_sum);
      row_data.PushBack( best_overlap_helix_sum);
      row_data.PushBack( best_overlap_strand_sum);
      row_data.PushBack( best_shift_sum);
      row_data.PushBack( best_shift_helix_sum);
      row_data.PushBack( best_shift_strand_sum);
      row_data.PushBack( best_length_dev_sum);
      row_data.PushBack( best_length_dev_helix_sum);
      row_data.PushBack( best_length_dev_strand_sum);
      row_data.PushBack( all_overlap_sum);
      row_data.PushBack( all_overlap_helix_sum);
      row_data.PushBack( all_overlap_strand_sum);
      row_data.PushBack( all_shift_sum);
      row_data.PushBack( all_shift_helix_sum);
      row_data.PushBack( all_shift_strand_sum);
      row_data.PushBack( all_length_dev_sum);
      row_data.PushBack( all_length_dev_helix_sum);
      row_data.PushBack( all_length_dev_strand_sum);

      // construct and insert row
      table.InsertRow( "pool", row_data);

      // end
      return table;
    }

    //! @brief initialize from a vector of sses
    //! @param SSE_VECTOR vector of sses
    //! @param IGNORE_UNSTRUCTURED ignore unstructured sses
    //! @param IDEALIZE whether to idealize the SSEs (or simply move them to the origin)
    void SSEPool::Initialize
    (
      const util::SiPtrVector< const SSE> &SSE_VECTOR,
      const bool IGNORE_UNSTRUCTURED,
      const bool IDEALIZE
    )
    {
      // iterate over sses in the vector
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( SSE_VECTOR.Begin()),
          sse_itr_end( SSE_VECTOR.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if not structured
        if( IGNORE_UNSTRUCTURED && !( *sse_itr)->GetType()->IsStructured())
        {
          continue;
        }

        // copy sse
        util::ShPtr< SSE> current_sse( ( *sse_itr)->Clone());

        // idealize if idealize param is true
        if( IDEALIZE)
        {
          // reset to unit conformation (ideal in origin)
          current_sse->SetToIdealConformationAtOrigin();
        }
        // just move to origin
        else
        {
          current_sse->Transform( math::Inverse( current_sse->GetOrientation()));
        }

        // insert into the set
        m_Data.Insert( current_sse);
      }
    }

  } // namespace assemble
} // namespace bcl
