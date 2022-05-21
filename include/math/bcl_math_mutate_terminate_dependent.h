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

#ifndef BCL_MATH_MUTATE_TERMINATE_DEPENDENT_H_
#define BCL_MATH_MUTATE_TERMINATE_DEPENDENT_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_interface.h"
#include "bcl_math_mutate_result.h"
#include "opti/bcl_opti_criterion_interface.h"
#include "opti/bcl_opti_tracker.h"
#include "random/bcl_random_distribution_interface.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateTerminateDependent
    //! @brief randomly selects from a list of mutates whose corresponding terminate criteria has not been met
    //! @details Contains a list of terminate criteria paired with a mutate object. Only mutate objects whose terminate
    //!          criteria has not been met are eligible to be used. Once the eligible mutate objects have been
    //!          determined, a random one is selected and used to mutate the t_Argument type.
    //!
    //! @see @link example_math_mutate_terminate_dependent.cpp @endlink
    //! @author alexanns
    //! @date Dec 13, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class MutateTerminateDependent :
      public MutateInterface< t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! list of terminates and paired mutate objects
      storage::List
      <
        storage::Pair
        <
          util::ShPtr< opti::CriterionInterface< t_ArgumentType, t_ResultType> >,
          util::ShPtr< MutateInterface< t_ArgumentType> >
        >
      > m_Mutates;

      //! scheme
      std::string m_Scheme;

      //! tracker
      util::SiPtr< const opti::Tracker< t_ArgumentType, t_ResultType> > m_Tracker;

      //! the random number generator used to select a random eligible mutate
      const random::DistributionInterface &m_RandomNumberGenerator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateTerminateDependent() :
        m_Mutates(),
        m_Scheme( GetStaticClassName< MutateTerminateDependent< t_ArgumentType, t_ResultType> >()),
        m_RandomNumberGenerator( random::GetGlobalRandom())
      {
      }

      //! @brief parameter constructor
      //! @param MUTATES list of terminates and paired mutate objects
      //! @param SCHEME scheme for this mutate object
      //! @param RNG the random number generator used to select a random eligible mutate
      MutateTerminateDependent
      (
        const storage::List
        <
          storage::Pair
          <
            util::ShPtr< opti::CriterionInterface< t_ArgumentType, t_ResultType> >,
            util::ShPtr< MutateInterface< t_ArgumentType> >
          >
        > &MUTATES,
        const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER,
        const std::string &SCHEME = ( GetStaticClassName< MutateTerminateDependent< t_ArgumentType, t_ResultType> >()),
        const random::DistributionInterface &RNG = random::GetGlobalRandom()
      ) :
        m_Mutates( MUTATES),
        m_Scheme( SCHEME),
        m_Tracker( TRACKER),
        m_RandomNumberGenerator( RNG)
      {
      }

      //! @brief Clone function
      //! @return pointer to new MutateTerminateDependent
      MutateTerminateDependent< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new MutateTerminateDependent< t_ArgumentType, t_ResultType>( *this);
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

      //! @brief return scheme
      //! @return scheme
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief set the tracker that this mutate points to
      //! @param TRACKER to tracker to reference for determining whether the criteria is met
      void SetTracker( const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER)
      {
        m_Tracker = util::ToSiPtr( TRACKER);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param ARGUMENT Argument of interest
      //! @return MutateResult that results from mutating to the argument
      MutateResult< t_ArgumentType> operator()( const t_ArgumentType &ARGUMENT) const
      {
        // true if there is nothing in the list of terminates and mutates
        if( m_Mutates.IsEmpty())
        {
          BCL_MessageStd( "There are no mutates in list");
          return MutateResult< t_ArgumentType>();
        }

        // vector to hold mutates that are eligible to be used - determined by their terminate criteria not being met
        storage::Vector< util::ShPtr< MutateInterface< t_ArgumentType> > > eligible_mutates;

        // build list of mutates that are eligible i.e. their terminate critieria has not been met yet
        for
        (
          typename storage::List
          <
            storage::Pair
            <
              util::ShPtr< opti::CriterionInterface< t_ArgumentType, t_ResultType> >,
              util::ShPtr< MutateInterface< t_ArgumentType> >
            >
          >::const_iterator itr( m_Mutates.Begin()), itr_end( m_Mutates.End());
          itr != itr_end; ++itr
        )
        {
          // true if the terminate has not had its criteria met yet - can add the mutate to the list of eligible mutates
          if( !itr->First()->CriteriaMet( *m_Tracker))
          {
            // add the mutate to the list of eligible mutates
            eligible_mutates.PushBack( itr->Second());
          }
        }

        // true if there is nothing in the list of mutates
        if( eligible_mutates.IsEmpty())
        {
          BCL_MessageStd( "No criteria were met");
          return MutateResult< t_ArgumentType>();
        }

        // get random index to mutate
        const size_t index( m_RandomNumberGenerator.Random< size_t>( 0, eligible_mutates.GetSize() - 1));

        // get the mutate corresponding to the index and use it to mutate "ARGUMENT"
        const MutateResult< t_ArgumentType> mutate_result( eligible_mutates( index)->operator()( ARGUMENT));

        // return the mutation result
        return mutate_result;
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
        // read members
        io::Serialize::Read( m_Mutates, ISTREAM);
        io::Serialize::Read( m_Scheme, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Mutates, OSTREAM, INDENT);
        io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // class MutateTerminateDependent

    // instantiate s_Instance
    template< typename t_DataType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> MutateTerminateDependent< t_DataType, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateTerminateDependent< t_DataType, t_ResultType>())
    );

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_MUTATE_TERMINATE_DEPENDENT_H_
