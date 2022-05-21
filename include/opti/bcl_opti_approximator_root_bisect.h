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

#ifndef BCL_OPTI_APPROXIMATOR_ROOT_BISECT_H_
#define BCL_OPTI_APPROXIMATOR_ROOT_BISECT_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_modular_base.h"
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorRootBisect
    //! @brief Implementation of the Bisection Method to find a root of a continuous function without a first
    //! derivative
    //! @details an objective function and an initial guess for a left and right side argument to the result are
    //! required
    //!
    //! @see http://en.wikipedia.org/wiki/Bisection_method
    //!
    //! @tparam t_ArgumentType an addition has to be defined for that type and a multiplication with double
    //! @tparam t_ResultType a '<' less than operator has to be defined
    //!
    //! @see @link example_opti_approximator_root_bisect.cpp @endlink
    //! @author woetzen, fischea
    //! @date Dec 18, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorRootBisect :
      public ApproximatorModularBase< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! ShPtr to the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > m_Objective;

      //! ShPtr to the left approximation result
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_Left;

      //! ShPtr to the right approximation result
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_Right;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorRootBisect() :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerAbsIsBetter),
        m_Objective(),
        m_Left(),
        m_Right()
      {
      }

      //! @brief construct from objective function and criterion
      //! @param OBJECTIVE the objective function
      //! @param LEFT_BORDER left border of the approximation interval
      //! @param RIGHT_BORDER right border of the approximation interval
      ApproximatorRootBisect
      (
        const math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> &OBJECTIVE,
        const CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const t_ArgumentType &LEFT_BORDER,
        const t_ArgumentType &RIGHT_BORDER
      ) :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerAbsIsBetter),
        m_Objective( util::CloneToShPtr( OBJECTIVE)),
        m_Left
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              LEFT_BORDER, m_Objective->operator()( LEFT_BORDER)
            )
          )
        ),
        m_Right
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              RIGHT_BORDER, m_Objective->operator()( RIGHT_BORDER)
            )
          )
        )
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorRootBisect< t_ArgumentType, t_ResultType>
      ApproximatorRootBisect *Clone() const
      {
        return new ApproximatorRootBisect( *this);
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

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief conducts the next approximation step and stores the approximation
      void Next()
      {
        // calculate next bisection
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_bisect( Bisect());

        // if left and middle function value have identical signs, replace left with middle
        if( ( m_Left->Second() < t_ResultType( 0)) == ( sp_bisect->Second() < t_ResultType( 0)))
        {
          m_Left = sp_bisect;
        }
        else // if right and middle function value have identical signs, replace right with middle
        {
          m_Right = sp_bisect;
        }

        // track current approximation
        this->Track( sp_bisect);
      }

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const
      {
        return true;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // call Read of base class
        ApproximatorModularBase< t_ArgumentType, t_ResultType>::Read( ISTREAM);

        // read members
        io::Serialize::Read( m_Objective, ISTREAM);
        io::Serialize::Read( m_Left, ISTREAM);
        io::Serialize::Read( m_Right, ISTREAM);

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
        io::Serialize::Write( m_Left, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Right, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief returns a ShPtr to the new argument and result pair in the middle
      //! @return ShPtr to the new argument and result pair in the middle
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > Bisect() const
      {
        // bisect the given interval
        const t_ArgumentType bisect( ( m_Left->First() + m_Right->First()) * 0.5);

        // return pair of argument and function value
        return util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
        (
          new storage::Pair< t_ArgumentType, t_ResultType>
          (
            bisect, m_Objective->operator()( bisect)
          )
        );
      }

    }; // template class ApproximatorRootBisect< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_ROOT_BISECT_H_
