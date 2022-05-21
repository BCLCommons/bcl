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

#ifndef BCL_SCORE_RESTRAINT_ATOM_ATTRACTION_H_
#define BCL_SCORE_RESTRAINT_ATOM_ATTRACTION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_trigonometric_transition.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintAtomAttraction
    //! @brief This class allows for an additional "Penalty" for atom distances that are too large or small by using a
    //!        piecewise function to put together a cos function with const functions at the top and bottom of the curve
    //!
    //! @see @link example_score_restraint_atom_attraction.cpp @endlink
    //! @author akinlr, weinerbe
    //! @date 06/30/2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintAtomAttraction :
      public RestraintAtomDistanceAssignment
    {
    private:

    //////////
    // data //
    //////////

      //! Well potential
      math::TrigonometricTransition m_WellPotential;

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

      //! @brief gives the default width of the cosine transition region of the peicewise function
      //! @return double which is the default width of the cosine transition region of the peicewise function
      static double GetDefaultTransitionWidth();

      //! @brief gives the default range of the depth of the cosine transition region of the peicewise function
      //! @return double which gives the default range of the depth of the cosine transition region of the function
      static const math::Range< double> &GetDefaultDepthRange();

      //! @brief amount to shift in the x-direction so that attraction is not level at same place as KB potential
      //! @return double amount to shift in x-direction so that attraction is not level at same place as KB potential
      static double GetDefaultScoreOffset();

      //! @brief gives the x-coordinate of the ending point for the attraction on the left of the KB potential
      //! @return double the x-coordinate of the ending point for the attraction on the left of the KB potential
      static double GetDefaultLeftEndWell( const math::Histogram &HISTOGRAM);

      //! @brief gives the x-coordinate of the ending point for the attraction on the right of the KB potential
      //! @return double the x-coordinate of the ending point for the attraction on the right of the KB potential
      static double GetDefaultRightEndWell( const math::Histogram &HISTOGRAM);

    //////////
    // data //
    //////////

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintAtomAttraction();

      //! @brief constructor from a specified histogram file
      //! @param WELL_DEPTH range over which the function will cover
      //! @param END_OF_WELL where the well should end
      //! @param WIDTH how wide the cos portion of the scoring function should be
      //! @param MAX_TO_MIN true if the function should go from max to min value - false otherwise
      //! @param SCHEME scheme to be used
      RestraintAtomAttraction
      (
        const math::Range< double> &WELL_DEPTH,
        const double END_OF_WELL,
        const double WIDTH,
        const bool MAX_TO_MIN,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone is the virtual copy constructor
      RestraintAtomAttraction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief gives the peicewise function that comprises this score
      //! @return function that is this score
      const math::TrigonometricTransition &GetFunction() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator scores protein model
      //! @param RESTRAINT restraint to be scored
      //! @return score
      double operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const;

    //////////////////////
    // input and output //
    //////////////////////

    private:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class RestraintAtomAttraction

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_ATOM_ATTRACTION_H_
