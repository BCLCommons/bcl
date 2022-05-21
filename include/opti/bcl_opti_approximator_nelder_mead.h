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

#ifndef BCL_OPTI_APPROXIMATOR_NELDER_MEAD_H_
#define BCL_OPTI_APPROXIMATOR_NELDER_MEAD_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_modular_base.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorNelderMead
    //! @brief ApproximatorNelderMead is an implementation of the Nelder-Mead method, a non linear optimization
    //! algorithm employing a simplex.
    //! @details It uses the concept of a simplex, which is a polytope of N + 1 vertices in N dimensions. It is a line
    //! segment on a line or a triangle on a plane and so on. It can find a locally optimal solution for a problem with
    //! N variables.
    //!
    //! @see http://en.wikipedia.org/wiki/Nelder-Mead_method
    //!
    //! @tparam t_ArgumentType the type of the argument
    //! @tparam t_ResultType the type of the result
    //!
    //! @see @link example_opti_approximator_nelder_mead.cpp @endlink
    //! @author fischea
    //! @date Dec 21, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorNelderMead :
      public ApproximatorModularBase< t_ArgumentType, t_ResultType>
    {

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class LessFunctionValue
      //! @brief Less struct to order list by t_ResultType
      //! @author fischea
      //! @remarks example unnecessary
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      struct LessFunctionValue
      {

        //! takes two arguments and returns true if A less B
        bool operator()
        (
          const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &RESULT_LHS,
          const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &RESULT_RHS
        ) const
        {
          return RESULT_LHS->Second() < RESULT_RHS->Second();
        }

      }; // struct LessFunctionValue

    //////////
    // data //
    //////////

      //! ShPtr to the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > m_Objective;

      //! ShPtrlist of arguments and results ordered by t_ResultType (small first)
      util::ShPtrList< storage::Pair< t_ArgumentType, t_ResultType> > m_Simplex;

      //! minimum difference between best and worst vertex in simplex
      t_ResultType m_TerminationDifference;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorNelderMead() :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerAbsIsBetter),
        m_Objective(),
        m_Simplex(),
        m_TerminationDifference()
      {
      }

      //! @brief construct from objective function and criterion
      //! @param OBJECTIVE the objective function
      //! @param CRITERION the termination criterion
      //! @param MINIMUM_DIFFERENCE minimum difference between best and worst vertex in simplex
      ApproximatorNelderMead
      (
        const math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> &OBJECTIVE,
        const CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const storage::List< t_ArgumentType> &START_SIMPLEX,
        const t_ResultType MINIMUM_DIFFERENCE
      ) :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerIsBetter),
        m_Objective( util::CloneToShPtr( OBJECTIVE)),
        m_TerminationDifference( MINIMUM_DIFFERENCE)
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);

        // fill the simplex
        FillSimplex( START_SIMPLEX);
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorNelderMead< t_ArgumentType, t_ResultType>
      ApproximatorNelderMead *Clone() const
      {
        return new ApproximatorNelderMead( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    private:

      //! @brief access to the best vertex in the simplex
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &GetBest() const
      {
        return m_Simplex.FirstElement();
      }

      //! @brief access to the second best vertex in the simplex
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &GetSecondBest() const
      {
        return *( ++m_Simplex.Begin());
      }

      //! @brief access to the worst vertex in the simplex
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &GetWorst() const
      {
        return m_Simplex.LastElement();
      }

      //! returns the reflection coefficient
      //! @return reflection coefficient
      static const double &GetReflectionCoefficient()
      {
        static const double s_reflection_coefficient( 1.0);
        return s_reflection_coefficient;
      }

      //! returns the expansion coefficient
      //! @return expansion coefficient
      static const double &GetExpansionCoefficient()
      {
        static const double s_expansion_coefficient( 2.0);
        return s_expansion_coefficient;
      }

      //! returns the outer contraction coefficient
      //! @return outer contraction coefficient
      static const double &GetOuterContractionCoefficient()
      {
        static const double s_outer_contraction_coefficient( 0.5);
        return s_outer_contraction_coefficient;
      }

      //! returns the inner contraction coefficient
      //! @return inner contraction coefficient
      static const double &GetInnerContractionCoefficient()
      {
        static const double s_inner_contraction_coefficient( -0.5);
        return s_inner_contraction_coefficient;
      }

      //! returns the shrink coefficient
      //! @return shrink coefficient
      static const double &GetShrinkCoefficient()
      {
        static const double s_shrink_coefficient( 0.5);
        return s_shrink_coefficient;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief conducts the next approximation step and stores the approximation
      void Next()
      {
        // references on best and second best (good)
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &best( GetBest());
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &good( GetSecondBest());

        // calculate center
        const t_ArgumentType center( Center());

        // reflect worst point along the center
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > reflected_value
        (
          Reflect( center, GetReflectionCoefficient())
        );

        // good point, which is better then the second best, but worse then the best value
        if( best->Second() <= reflected_value->Second() && reflected_value->Second() <= good->Second())
        {
          //replace the worst with the reflected value
          Replace( reflected_value);
        }
        // if reflection turns out better than the best, try to make a larger reflection
        else if( reflected_value->Second() < best->Second())
        {
          // try to expand reflection
          const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > expand_value
          (
            Reflect( center, GetExpansionCoefficient())
          );

          // if expansion is still better the the best value, take the expanded, but try not to go further
          if( expand_value->Second() < reflected_value->Second())
          {
            //replace worst with expansion
            Replace( expand_value);
          }
          // expansion did not help, just take the initial best value
          else
          {
            // replace the worst with the reflected value
            Replace( reflected_value);
          }
        }
        else // reflection did not improve
        {
          // try outer reflection
          const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > outer_value
          (
            Reflect( center, GetOuterContractionCoefficient())
          );

          // outer reflection better then second best value
          if( outer_value->Second() < good->Second())
          {
            Replace( outer_value);
          }
          else
          {
            // try inner reflection
            const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > inner_value
            (
              Reflect( center, GetInnerContractionCoefficient())
            );

            // inner reflection better then second best value
            if( inner_value->Second() < good->Second())
            {
              Replace( inner_value);
            }
            else
            {
              // shrink complete simplex
              Shrink();
            }
          }
        }

        this->Track( GetBest());
      }

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const
      {
        return !( GetWorst()->Second() - GetBest()->Second() < m_TerminationDifference);
      }

    private:

      //! takes a storage list of arguments and calculates the values to all arguments and puts them in the list
      //! @param SIMPLEX argument list for the simplex to fill
      void FillSimplex( const storage::List< t_ArgumentType> &SIMPLEX)
      {
        // reset the simplex
        m_Simplex.Reset();

        // iterate over all arguments
        for
        (
          typename storage::List< t_ArgumentType>::const_iterator arg_itr( SIMPLEX.Begin());
          arg_itr != SIMPLEX.End();
          ++arg_itr
        )
        {
          // push back a pair of the current argument and the function values for this argument in the list
          m_Simplex.PushBack
          (
            util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
            (
              new storage::Pair< t_ArgumentType, t_ResultType>
              (
                *arg_itr,
                m_Objective->operator()( *arg_itr)
              )
            )
          );
        }

        // sort the simplex by function value
        m_Simplex.Sort( LessFunctionValue());
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
        // call Read of base class
        ApproximatorModularBase< t_ArgumentType, t_ResultType>::Read( ISTREAM);

        // read members
        io::Serialize::Read( m_Objective, ISTREAM);
        io::Serialize::Read( m_Simplex, ISTREAM);
        io::Serialize::Read( m_TerminationDifference, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // call Write of base class
        ApproximatorModularBase< t_ArgumentType, t_ResultType>::Write( OSTREAM, INDENT) << '\n';

        // write members
        io::Serialize::Write( m_Objective, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Simplex, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_TerminationDifference, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! replaces a pair of argument and result by a better one
      //! @param PAIR_INSERT ShPtr to the pair to be inserted
      void Replace( const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &PAIR_INSERT)
      {
        // remove the worst from the end of the simplex
        m_Simplex.PopBack();

        // insert the new one (where does not matter, since the simplex has to be resorted
        m_Simplex.PushBack( PAIR_INSERT);

        // sort the simplex by the Result, smallest first
        m_Simplex.Sort( LessFunctionValue());
      }

      //! calculate center of points without worst point
      //! @return center of points without worst point
      t_ArgumentType Center() const
      {
        // copy the best argument as the center
        t_ArgumentType center( GetBest()->First());

        // iterate over all points worse then the best one excluding the worst, which will be reflected
        for
        (
          typename util::ShPtrList< storage::Pair< t_ArgumentType, t_ResultType> >::const_iterator
            pair_itr( ++m_Simplex.Begin()), pair_itr_end( --m_Simplex.End());
          pair_itr != pair_itr_end;
          ++pair_itr
        )
        {
          // add point to the best argument
          center += ( *pair_itr)->First();
        }

        // divide by number of points minus the worst, which will be reflected
        center /= double( m_Simplex.GetSize() - 1);

        return center;
      }

      //! @brief returns ShPtr to the reflected worst point
      //! @param CENTER the worst point and the center define a line, on which the reflection is happening
      //! @param COEFFICIENT determines the direction and magnitude of the reflection
      //! @return ShPtr to the reflected worst point
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
      Reflect( const t_ArgumentType &CENTER, const double &COEFFICIENT) const
      {
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > reflected_worst_point
        (
          new storage::Pair< t_ArgumentType, t_ResultType>
          (
            CENTER + COEFFICIENT * ( CENTER - m_Simplex.LastElement()->First()),
            t_ResultType()
          )
        );

        reflected_worst_point->Second() = m_Objective->operator()( reflected_worst_point->First());

        return reflected_worst_point;
      }

      //! shrink every distance to best point by m_Shrink_coefficient
      void Shrink()
      {
        // make a reference to best point
        const t_ArgumentType &best_point( GetBest()->First());
        util::ShPtrList< storage::Pair< t_ArgumentType, t_ResultType> > shrunk_list;
        shrunk_list.PushFront( GetBest());

        // iterate over all points worse then the best one
        for
        (
          typename util::ShPtrList< storage::Pair< t_ArgumentType, t_ResultType> >::const_iterator
            pair_itr( ++m_Simplex.Begin()),
            pair_itr_end( m_Simplex.End());
          pair_itr != pair_itr_end;
          ++pair_itr
        )
        {
          // calculate shrunken point
          const t_ArgumentType new_point( best_point + GetShrinkCoefficient() * ( ( *pair_itr)->First() - best_point));

          // store shrunken point in new list
          shrunk_list.PushBack
          (
            util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
            (
              new storage::Pair< t_ArgumentType, t_ResultType>
              (
                new_point, m_Objective->operator()( new_point)
              )
            )
          );
        }

        // replace old simplex with shrunken simplex
        m_Simplex = shrunk_list;
        m_Simplex.Sort( LessFunctionValue());
      }

    }; // template class ApproximatorNelderMead< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_NELDER_MEAD_H_
