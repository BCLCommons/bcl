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

#ifndef BCL_OPTI_APPROXIMATOR_GOLDEN_SECTION_H_
#define BCL_OPTI_APPROXIMATOR_GOLDEN_SECTION_H_

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
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorGoldenSection
    //! @brief Implementation of the golden section search method to find a minimum of a unimodal function
    //! @details and objective function and two initial arguments are provided, that are limiting the solution
    //! the iterate might return the same best result for multiple iterations, but the interval will shrink
    //! while the confidence in the result will increase
    //!
    //! @see http://en.wikipedia.org/wiki/Golden_section_search
    //!
    //! @tparam t_ArgumentType an addition and subtraction, as well as a scalar multiplication has to be defined
    //! @tparam t_ResultType a less than comparison has to be defined
    //!
    //! @see @link example_opti_approximator_golden_section.cpp @endlink
    //! @author woetzen, fischea
    //! @date Dec 18, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorGoldenSection :
       public ApproximatorModularBase< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! ShPtr to the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > m_Objective;

      //! left argument result pair
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_Left;

      //! right argument result pair
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_Right;

      //! middle argument result pair
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_Middle;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorGoldenSection() :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerIsBetter),
        m_Objective(),
        m_Left(),
        m_Right(),
        m_Middle()
      {
      }

      //! @brief construct from objective function and criterion
      //! @param OBJECTIVE the objective function
      //! @param CRITERION the termination criterion
      //! @param ARGUMENT initial argument
      ApproximatorGoldenSection
      (
        const math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> &OBJECTIVE,
        const CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const t_ArgumentType &BORDER_LEFT,
        const t_ArgumentType &BORDER_RIGHT
      ) :
        ApproximatorModularBase< t_ArgumentType, t_ResultType>( e_SmallerIsBetter),
        m_Objective( util::CloneToShPtr( OBJECTIVE)),
        m_Left
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              BORDER_LEFT, m_Objective->operator()( BORDER_LEFT)
            )
          )
        ),
        m_Right
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              BORDER_RIGHT, m_Objective->operator()( BORDER_RIGHT)
            )
          )
        ),
        m_Middle( GoldenIntersect( BORDER_LEFT, BORDER_RIGHT))
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);

        // track the initial approximation
        this->Track( m_Left);
        this->Track( m_Middle);
        this->Track( m_Right);
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorGoldenSection< t_ArgumentType, t_ResultType>
      ApproximatorGoldenSection *Clone() const
      {
        return new ApproximatorGoldenSection( *this);
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

      //! @brief returns 2 - golden ratio
      //! @return 2 - golden ratio
      static const double &GetResPhi()
      {
        static const double s_res_phi( 2.0 - ( 1.0 + math::Sqrt( 5.0)) * 0.5);
        return s_res_phi;
      }

    protected:

      //! @brief conducts the next approximation step and stores the approximation
      void Next()
      {
        // device interval for next pair of argument and function value
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_next
        (
          GoldenIntersect( m_Middle->First(), m_Right->First())
        );

        // get best approximation
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_best
        (
          this->GetTracker().GetBest()
        );

        // update other, if the next is better
        if( sp_next->Second() < sp_best->Second())
        {
          m_Left = m_Middle;
          m_Middle = sp_next;
          // m_Right stays the same
        }
        else
        {
          // m_Middle stays the same
          m_Right = m_Left; // switch right and left
          m_Left = sp_next;
        }

        // update best
        this->Track( sp_next);
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
        io::Serialize::Read( m_Left, ISTREAM);
        io::Serialize::Read( m_Right, ISTREAM);
        io::Serialize::Read( m_Middle, ISTREAM);

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
        io::Serialize::Write( m_Right, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Middle, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief intersects the given interval with the golden ratio and returns the objective function
      //! value at this point
      //! @param LEFT left border of the interval
      //! @param RIGHT right border of the interval
      //! @return objective function value at the intersection point
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > GoldenIntersect
      (
        const t_ArgumentType &LEFT,
        const t_ArgumentType &RIGHT
      ) const
      {
        // device interval for next pair of argument and function value
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_intersect
        (
          new storage::Pair< t_ArgumentType, t_ResultType>
          (
            LEFT + GetResPhi() * ( RIGHT - LEFT), t_ResultType()
          )
        );
        sp_intersect->Second() = m_Objective->operator()( sp_intersect->First());

        return sp_intersect;
      }

    }; // template class ApproximatorGoldenSection< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_GOLDEN_SECTION_H_
