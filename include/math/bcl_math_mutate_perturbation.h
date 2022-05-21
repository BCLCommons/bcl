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

#ifndef BCL_MATH_MUTATE_PERTURBATION_H_
#define BCL_MATH_MUTATE_PERTURBATION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_interface.h"
#include "bcl_math_mutate_result.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutatePerturbation
    //! @brief template class MutatePerturbation
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_math_mutate_perturbation.cpp @endlink
    //! @author woetzen
    //! @date Aug 23, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class MutatePerturbation :
      public util::ObjectInterface
    {

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class SinglePerturbation
      //! @brief a single perturbation in the mutate
      //! @remarks example unnecessary
      //! @author woetzen
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class SinglePerturbation :
        public util::ObjectInterface
      {
        friend class MutatePerturbation;
      private:

      //////////
      // data //
      //////////

        size_t                        m_Total_repetition;
        size_t                        m_Current_repetition;
        t_ArgumentType                m_CurrentArgument;
        util::ShPtr< MutateInterface< t_ArgumentType> > m_Perturbation;

      //////////////////////////////////
      // construction and destruction //
      //////////////////////////////////

      public:

        //! @brief defautl constructor
        SinglePerturbation() :
          m_Total_repetition( 0),
          m_Current_repetition( 0),
          m_CurrentArgument(),
          m_Perturbation()
        {
        }

        //! @brief construct from Perturbation and total repetitions and argument
        //! @param SP_MUTATE ShPtr to mutate defining the perturbation
        //! @param REPETITION the total number of repetitions that need to be performed
        //! @param ARGUMENT the starting argument
        SinglePerturbation
        (
          const util::ShPtr< MutateInterface< t_ArgumentType> > &SP_MUTATE,
          const size_t REPETITION,
          const t_ArgumentType &ARGUMENT
        ) :
          m_Total_repetition( REPETITION),
          m_Current_repetition( 0),
          m_CurrentArgument( ARGUMENT),
          m_Perturbation( SP_MUTATE)
        {
        }

        //! @brief Clone function
        //! @return pointer to new SinglePerturbation
        SinglePerturbation *Clone() const
        {
          return new SinglePerturbation( *this);
        }

        //! @brief returns class name of the object behind a pointer or the current object
        //! @return the class name
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

      private:

      ////////////////
      // operations //
      ////////////////

        //! @brief reset perturbation
        void Reset( const t_ArgumentType &ARGUMENT)
        {
          m_Current_repetition = 0;
          m_CurrentArgument = ARGUMENT;
        }

        //! @brief perturb current argument
        void Perturb()
        {
          // store the result
          MutateResult< t_ArgumentType> mutate_result( m_Perturbation->operator()( m_CurrentArgument));

          // if the argument is not defined
          while( !mutate_result.GetArgument().IsDefined())
          {
            // do not change the current argument but issue a warning
            mutate_result = m_Perturbation->operator()( m_CurrentArgument);
          }
          m_CurrentArgument = *mutate_result.GetArgument();

          // increment the repetition number
          ++m_Current_repetition;
        }

        //! @brief is perturbation finished?
        //! @return true if current repetition meets or is higher than total repetitions
        bool IsFinished() const
        {
          return m_Current_repetition >= m_Total_repetition;
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
          // read member
          io::Serialize::Read( m_Total_repetition, ISTREAM);
          io::Serialize::Read( m_Current_repetition, ISTREAM);
          io::Serialize::Read( m_CurrentArgument, ISTREAM);
          io::Serialize::Read( m_Perturbation, ISTREAM);

          // end
          return ISTREAM;
        }

        //! @brief write to std::ostream
        //! @param OSTREAM output stream to write to
        //! @param INDENT indentation
        //! @return output stream which was written to
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
        {
          // write member
          io::Serialize::Write( m_Total_repetition, OSTREAM, INDENT) << '\n';
          io::Serialize::Write( m_Current_repetition, OSTREAM, INDENT) << '\n';
          io::Serialize::Write( m_CurrentArgument, OSTREAM, INDENT) << '\n';
          io::Serialize::Write( m_Perturbation, OSTREAM, INDENT);

          // end
          return OSTREAM;
        }

      }; // class SinglePerturbation

    private:

    //////////
    // data //
    //////////

      t_ArgumentType                                          m_StartArgument;
      storage::Vector< SinglePerturbation>                    m_Perturbations;
      typename storage::Vector< SinglePerturbation>::iterator m_CurrentPerturbation;

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
      MutatePerturbation< t_ArgumentType>()
      {
      }

      //! @brief construct from starting argument
      MutatePerturbation< t_ArgumentType>( const t_ArgumentType &START) :
        m_StartArgument( START),
        m_Perturbations( 0),
        m_CurrentPerturbation( m_Perturbations.Begin())
      {
      }

      //! @brief Clone function
      //! @return pointer to new MutatePerturbation< t_ArgumentType>
      MutatePerturbation< t_ArgumentType> *Clone() const
      {
        return new MutatePerturbation< t_ArgumentType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get start argument
      const t_ArgumentType &GetStartArgument() const
      {
        return m_StartArgument;
      }

      //! @brief set the start argument
      //! resets perturbation and replaces the start argument with the given one
      //! @param ARGUMENT start argument to be used
      void SetStartArgument( const t_ArgumentType &ARGUMENT)
      {
        m_StartArgument = ARGUMENT;
        Reset();
      }

      //! @brief get current argument
      //! @return current argument
      const t_ArgumentType &GetCurrent() const
      {
        return m_CurrentPerturbation->m_CurrentArgument;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add an additionla perturbation level
      //! @param MUTATE mutation to be performed
      //! @param REPETITION how often to repeat
      void InsertPerturbation( const util::ShPtr< MutateInterface< t_ArgumentType> > &MUTATE, const size_t REPETITION)
      {
        m_Perturbations.PushBack( SinglePerturbation( MUTATE, REPETITION, m_StartArgument));
        Reset();
      }

      //! @brief reset perturbation
      void Reset()
      {
        // reset the current perturbation to the start
        m_CurrentPerturbation = m_Perturbations.Begin();

        // reset all single perturbations
        for
        (
          typename storage::Vector< SinglePerturbation>::iterator itr( m_Perturbations.Begin()), itr_end( m_Perturbations.End());
          itr != itr_end;
          ++itr
        )
        {
          itr->Reset( m_StartArgument);
        }
      }

      //! @brief check if perturbations are finished
      //! @return true is last iteration is finished
      bool IsFinished() const
      {
        return m_Perturbations.LastElement().IsFinished();
      }

      //! @brief perturb
      //! @return the current perturbation
      const t_ArgumentType &Perturb()
      {
        // all iterations are finished if last element is finished
        if( IsFinished())
        {
          Reset();
        }

        if( m_CurrentPerturbation->IsFinished())
        {
          // go to first not finished iteration
          do
          {
            ++m_CurrentPerturbation;
          }
          while( m_CurrentPerturbation->IsFinished());

          // perturb
          m_CurrentPerturbation->Perturb();

          // set all previously finished iteration to the result of that iteration
          for
          (
            typename storage::Vector< SinglePerturbation>::iterator itr( m_Perturbations.Begin());
            itr != m_CurrentPerturbation;
            ++itr
          )
          {
            itr->Reset( m_CurrentPerturbation->m_CurrentArgument);
          }
          // set the current perturbation to the start
          m_CurrentPerturbation = m_Perturbations.Begin();
        }
        else
        {
          // perturb
          m_CurrentPerturbation->Perturb();
        }

        // end
        return GetCurrent();
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
        // read member
        io::Serialize::Read( m_StartArgument, ISTREAM);
        io::Serialize::Read( m_Perturbations, ISTREAM);

        // read current perturbation
        size_t current( 0);
        io::Serialize::Read( current, ISTREAM);
        m_CurrentPerturbation = m_Perturbations.Begin() + current;

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_StartArgument, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Perturbations, OSTREAM, INDENT) << '\n';

        // write current perturbation
        const size_t current_perturbation( m_CurrentPerturbation - m_Perturbations.Begin());
        io::Serialize::Write( current_perturbation, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class MutatePerturbation

    // instantiate s_Instance
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> MutatePerturbation< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new MutatePerturbation< t_ArgumentType>())
    );

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_MUTATE_PERTURBATION_H_
