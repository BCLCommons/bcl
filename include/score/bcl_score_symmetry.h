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

#ifndef BCL_SCORE_SYMMETRY_H_
#define BCL_SCORE_SYMMETRY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "coord/bcl_coord_movable_interface.h"
#include "find/bcl_find_locator_interface.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Symmetry
    //! @brief scores the symmetry of an object.
    //!
    //! @see @link example_score_symmetry.cpp @endlink
    //! @author alexanns
    //! @date June 26, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class Symmetry :
      public math::FunctionInterfaceSerializable< t_ArgumentType, double>
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used in outputting schemes
      std::string m_Scheme;

      storage::Vector
      <
        util::ShPtrVector< find::LocatorInterface< linal::Vector3D, t_ArgumentType> >
      > m_MovableLocators;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme()
      {
        // static string
        static const std::string s_default_scheme( "sym");

        // end
        return s_default_scheme;
      }

      //! @brief return command line flag for using symmetry score
      //! @return command line flag for using symmetry score
      static util::ShPtr< command::FlagInterface> &GetFlagScoreSymmetry()
      {
        // initialize static flag
        static util::ShPtr< command::FlagInterface> s_flag
        (
          new command::FlagStatic
          (
            "score_symmetry", "\tflag to enable use of symmetry scores, requires a symmetry definition file and a weight",
            command::Parameter
            (
              "symmetry_definition_filename",
              "\tfull path and name of file containing the symmetry definitions",
              ""
            )
          )
        );
        // end
        return s_flag;
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking member argument
      Symmetry
      (
        const storage::Vector
        <
          util::ShPtrVector< find::LocatorInterface< linal::Vector3D, t_ArgumentType> >
        > &MOVABLE_LOCATORS,
        const std::string &SCHEME = GetDefaultScheme()
      ) :
        m_Scheme( SCHEME),
        m_MovableLocators( MOVABLE_LOCATORS)
      {
      }

      //! @brief construct from file that contains movable locators
      //! this function has to be specialized for each t_ArgumentType within the cpp
      //! @param FILENAME filename that contains some kind of movable locators
      //! @param SCHEME
      Symmetry
      (
        const std::string &FILENAME,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new Counter
      Symmetry< t_ArgumentType> *Clone() const
      {
        return new Symmetry< t_ArgumentType>( *this);
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

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() used for scoring the symmetry
      //! @param ARGUMENT the t_Argument which the will be scored for being symmetric
      double operator()( const t_ArgumentType &ARGUMENT) const
      {
        // create "coordinates" to hold the coordinates of all the symmetric units
        util::ShPtrVector< util::SiPtrVector< const linal::Vector3D> > coordinates( m_MovableLocators.GetSize(), util::SiPtrVector< const linal::Vector3D>());
        util::ShPtrVector< linal::Vector3D> all_coordinates;

        size_t num_symmetric_units( m_MovableLocators.GetSize());
        for( size_t position( 0); position < m_MovableLocators.Begin()->GetSize(); ++position)
        {
          util::SiPtrVector< const linal::Vector3D> current_symmetric_units_points;
          bool symmetry_point( true);

          // get the symmetry point for each symmetric unit
          for
          (
            size_t symmetric_unit( 0);
            symmetric_unit < num_symmetric_units && symmetry_point;
            ++symmetric_unit
          )
          {
            const linal::Vector3D movable( m_MovableLocators( symmetric_unit)( position)->Locate( ARGUMENT));

            // check if this movable interface actually exists
            if( !movable.IsDefined())
            {
              symmetry_point = false;
              break;
            }
            all_coordinates.PushBack( util::ShPtr< linal::Vector3D>( movable.Clone()));
            const util::SiPtr< const linal::Vector3D> point( *all_coordinates.LastElement());

            BCL_MessageDbg( "vector 3d is " + util::Format()( *point));

            //
            symmetry_point &= point->IsDefined();

            // add the current symmetry point coordinates for the unit denoted by "itr" to "current_symmetric_units_points"
            current_symmetric_units_points.PushBack( point);
          }

          // all of the current symmetry point coordinates were defined so can add all of them to "coordinates"
          if( symmetry_point)
          {
            // iterate through "current_symmetric_units_points" to add it to the appropriate siptrvector in "coordinates"
            for
            (
              util::SiPtrVector< const linal::Vector3D>::const_iterator
                sym_unit_itr( current_symmetric_units_points.Begin()),
                sym_unit_iter_end( current_symmetric_units_points.End());
              sym_unit_itr != sym_unit_iter_end;
              ++sym_unit_itr
            )
            {
              coordinates( sym_unit_itr - current_symmetric_units_points.Begin())->PushBack( *sym_unit_itr);
            }
          }
//          else
//          {
//            BCL_MessageDbg( "found points for symmetry layer, but at least one is undefined");
//          }
        }

        double score( coord::SymmetryFactor( coordinates));
//        if( score != 0)
//        {
//          score = -1.0/score;
//        }

        if( !util::IsDefined( score))
        {
          score = 0.0;
        }

        BCL_MessageStd
        (
          "symmetry score considers " + util::Format()( coordinates.GetSize()) + " sets of coordinates with " +
          util::Format()( coordinates( 0)->GetSize()) + " elements each."
        );
        BCL_MessageStd( "Final symmetry score is " + util::Format()( score));

        return score;
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
        io::Serialize::Read( m_MovableLocators, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_MovableLocators, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // class Symmetry

    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> Symmetry< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Symmetry< t_ArgumentType>())
    );

    //! @brief construct from file that contains movable locators
    //! this function has to be specialized for each t_ArgumentType within the cpp
    //! @param FILENAME filename that contains some kind of movable locators
    //! @param SCHEME
    template<>
    BCL_API Symmetry< assemble::ProteinModel>::Symmetry
    (
      const std::string &FILENAME,
      const std::string &SCHEME
    );

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_SYMMETRY_H_ 
