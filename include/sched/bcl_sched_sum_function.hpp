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

#ifndef BCL_SCHED_SUM_FUNCTION_HPP_
#define BCL_SCHED_SUM_FUNCTION_HPP_

// include header of this class
#include "bcl_sched_sum_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_sched_function_job_with_data.h"
#include "bcl_sched_scheduler_interface.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sched
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> SumFunction< t_ArgumentType, t_ResultType>::s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_ArgumentType, typename t_ResultType>
    SumFunction< t_ArgumentType, t_ResultType>::SumFunction( const std::string &SCHEME) :
      Base( SCHEME)
    {
    }

    //! @brief construct from ShPtr on Function, optional coefficient, and absolute value
    //! y = y0 + s1*g1(x)
    //! @param ABSOLUTE    y0 as absolute in sumfunction - default 0
    //! @param COEFFICIENT s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
    //! @param SP_FUNCTION g1(x) as part of sumfunction as ShPtr to FunctionInterface
    template< typename t_ArgumentType, typename t_ResultType>
    SumFunction< t_ArgumentType, t_ResultType>::SumFunction
    (
      const t_ResultType &ABSOLUTE,
      const double &COEFFICIENT,
      const util::Implementation< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > &SP_FUNCTION,
      const std::string &SCHEME
    ) :
      Base( ABSOLUTE, COEFFICIENT, SP_FUNCTION, SCHEME)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new SumFunction which is a copy of this
    template< typename t_ArgumentType, typename t_ResultType>
    SumFunction< t_ArgumentType, t_ResultType> *SumFunction< t_ArgumentType, t_ResultType>::Clone() const
    {
      return new SumFunction< t_ArgumentType, t_ResultType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_ArgumentType, typename t_ResultType>
    const std::string &SumFunction< t_ArgumentType, t_ResultType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief actual f(x) operator taking ARGUMENT and returning the sum of all function values
    //! @brief ARGUMENT x to f(x)
    //! @return the value of f(x) = sum
    template< typename t_ArgumentType, typename t_ResultType>
    t_ResultType
    SumFunction< t_ArgumentType, t_ResultType>::operator()
    (
      const t_ArgumentType &ARGUMENT
    ) const
    {
      // list to store jobs an their result
      storage::List< storage::Triplet< util::ShPtr< JobInterface>, t_ResultType, double> > jobs_result;

      // iterate over all terms
      for( const_iterator itr( this->GetFunction().Begin()), itr_end( this->GetFunction().End()); itr != itr_end; ++itr)
      {
        // create new job trans_center result triplet
        jobs_result.PushBack( storage::Triplet< util::ShPtr< JobInterface>, t_ResultType, double>());

        // create reference to the just create triplet
        storage::Triplet< util::ShPtr< JobInterface>, t_ResultType, double> &current_job_result( jobs_result.LastElement());

        // initialize the actual job object
        current_job_result.First() =
        util::ShPtr< JobInterface>
        (
          new FunctionJobWithData< t_ArgumentType, t_ResultType>
          (
            0,
            util::ToSiPtr( *itr->Second()),
            ARGUMENT,
            JobInterface::e_READY,
            &current_job_result.Second()
          )
        );
        current_job_result.Third() = itr->First(); // add the scalar

        // submit job to the scheduler
        GetScheduler().SubmitJob( current_job_result.First());
      }

      t_ResultType sum( this->GetAbsolute());
      // iterate over jobs, join and add results
      for
      (
        typename storage::List< storage::Triplet< util::ShPtr< JobInterface>, t_ResultType, double> >::iterator
          itr( jobs_result.Begin()), itr_end( jobs_result.End());
        itr != itr_end;
        ++itr
      )
      {
        GetScheduler().Join( itr->First());
        sum += itr->Third() * itr->Second();
      }

      // end
      return sum;
    }

  } // namespace sched
} // namespace bcl

#endif // BCL_SCHED_SUM_FUNCTION_HPP_
