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
#include "density/bcl_density_simulate_default.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SimulateDefault::s_Instance
    (
      GetObjectInstances().AddInstance( new SimulateDefault( linal::Vector3D( 2.3, 2.3, 2.3), 6.9, e_Gaussian))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from grid spacing and resolution
    //! @param GRID_SPACING the spacing for the density grid
    //! @param RESOLUTION the resolution to simulate for
    //! @param SMOOTHING_KERNEL kernel for the smoothing
    SimulateDefault::SimulateDefault
    (
      const linal::Vector3D &GRID_SPACING,
      const double RESOLUTION,
      const Kernel SMOOTHING_KERNEL
    ) :
      m_GridSpacing( GRID_SPACING),
      m_Resolution( RESOLUTION),
      m_Margin( 2),
      m_SmoothingKernel( SMOOTHING_KERNEL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SimulateDefault
    SimulateDefault *SimulateDefault::Clone() const
    {
      return new SimulateDefault( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SimulateDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution
    //! @param RESOLUTION the resolution for the density map to be generated
    void SimulateDefault::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set the resolution
    double SimulateDefault::GetResolution() const
    {
      return m_Resolution;
    }

    //! @brief set the grid spacing
    //! @param GRID_SPACING the width of a grid element in x, y and z
    void SimulateDefault::SetGridSpacing( const linal::Vector3D &GRID_SPACING)
    {
      m_GridSpacing = GRID_SPACING;
    }

    //! @brief set the margin
    //! @param MARGIN number of additional cells next to last atom occupied cells
    void SimulateDefault::SetMargin( const size_t MARGIN)
    {
      m_Margin = MARGIN;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief generate simulated density from given list of atoms
    //! @param ATOMS siptrvector of atoms
    //! @return a simulated density map
    Map SimulateDefault::operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const
    {
      // check that any atoms are given
      if( ATOMS.IsEmpty())
      {
        return Map();
      }

      // measure ATOMS extent
      storage::VectorND< 2, linal::Vector3D> min_max_coord( DetermineGridCorners( ATOMS));

      linal::Vector3D &mincoord( min_max_coord.First());
      linal::Vector3D &maxcoord( min_max_coord.Second());

      // bring lattice into register with origin
      for( size_t i( 0); i < 3; ++i)
      {
        mincoord( i) = m_GridSpacing( i) * std::floor( mincoord( i) / m_GridSpacing( i));
        maxcoord( i) = m_GridSpacing( i) * std::ceil( maxcoord( i) / m_GridSpacing( i));
      }

      // allocate atom density map
      const storage::VectorND< 3, size_t> ext_atom_map // extension of atom map
      (
        size_t( ceil( ( maxcoord.X() - mincoord.X()) / m_GridSpacing.X()) + 2 * m_Margin + 1),
        size_t( ceil( ( maxcoord.Y() - mincoord.Y()) / m_GridSpacing.Y()) + 2 * m_Margin + 1),
        size_t( ceil( ( maxcoord.Z() - mincoord.Z()) / m_GridSpacing.Z()) + 2 * m_Margin + 1)
      );
      const linal::Vector3D grid2
      (
        mincoord.X() - m_Margin * m_GridSpacing.X(),
        mincoord.Y() - m_Margin * m_GridSpacing.Y(),
        mincoord.Z() - m_Margin * m_GridSpacing.Z()
      );

      math::Tensor< double> atom_map( ext_atom_map.Third(), ext_atom_map.Second(), ext_atom_map.First(), double( 0));

      // correct for lattice interpolation smoothing effects
      // slightly lowers the kernel width to maintain target resolution
      const size_t corrmode( 2);

      // desired kernel amplitude (scaling factor):
      const double kernel_amplitude( 1.0);

      // interpolate structure to protein map and keep track of variability
      // Projecting atoms to cubic lattice by trilinear interpolation...
      double varp( 0.0);
      size_t total_weight( 0);
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator
          atom_itr( ATOMS.Begin()), atom_itr_end( ATOMS.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        const biol::Atom &current_atom( **atom_itr);
        const linal::Vector3D &current_coord( current_atom.GetCoordinates());

        if( !current_coord.IsDefined())
        {
          continue;
        }
        // weigh by atom type
        size_t weight( 0);
        if( !current_atom.GetType()->GetElementType().IsDefined())
        {
          continue;
        }
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Lithium)  continue;   // till Lithium
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Neon)     weight = 1; // till Neon
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Sulfur)   weight = 2; // till Sulfur
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Titanium) weight = 3; // till Titanium
        else if( current_atom.GetType()->GetElementType() <= chemistry::GetElementTypes().e_Nickel)   weight = 4; // till Nickel
        else weight = 5;                                             // from Copper

        // sum up weights
        total_weight += weight;

        // compute position within grid in Angstroem
        const linal::Vector3D grid_position
                             (
                               m_Margin + ( current_coord.X() - mincoord.X()) / m_GridSpacing.X(),
                               m_Margin + ( current_coord.Y() - mincoord.Y()) / m_GridSpacing.Y(),
                               m_Margin + ( current_coord.Z() - mincoord.Z()) / m_GridSpacing.Z()
                             );
        // convert to
        const int x0( int( std::floor( grid_position.X())));
        const int y0( int( std::floor( grid_position.Y())));
        const int z0( int( std::floor( grid_position.Z())));

        // interpolate - determine position in the box that point is interpolated to
        const double a( x0 + 1 - grid_position.X());
        const double b( y0 + 1 - grid_position.Y());
        const double c( z0 + 1 - grid_position.Z());

        const double ab( a * b);
        const double ac( a * c);
        const double bc( b * c);
        const double aa( a * a);
        const double bb( b * b);
        const double cc( c * c);
        const double oma( 1 - a);
        const double omb( 1 - b);
        const double omc( 1 - c);
        const double omab( oma * omb);
        const double omac( oma * omc);
        const double ombc( omb * omc);
        const double omaa( oma * oma);
        const double ombb( omb * omb);
        const double omcc( omc * omc);
        const double val1(   ab *    c * weight);
        const double val2(   ab *  omc * weight);
        const double val3(   ac *  omb * weight);
        const double val4(   bc *  oma * weight);
        const double val5(    a * ombc * weight);
        const double val6(    c * omab * weight);
        const double val7(    b * omac * weight);
        const double val8( omab * omc  * weight);

        atom_map( z0    , y0    , x0    ) += val1;
        atom_map( z0 + 1, y0    , x0    ) += val2;
        atom_map( z0    , y0 + 1, x0    ) += val3;
        atom_map( z0    , y0    , x0 + 1) += val4;
        atom_map( z0 + 1, y0 + 1, x0    ) += val5;
        atom_map( z0    , y0 + 1, x0 + 1) += val6;
        atom_map( z0 + 1, y0    , x0 + 1) += val7;
        atom_map( z0 + 1, y0 + 1, x0 + 1) += val8;

        varp += val1 * ( omaa + ombb + omcc);
        varp += val2 * ( omaa + ombb +   cc);
        varp += val3 * ( omaa +   bb + omcc);
        varp += val4 * (   aa + ombb + omcc);
        varp += val5 * ( omaa +   bb +   cc);
        varp += val6 * (   aa +   bb + omcc);
        varp += val7 * (   aa + ombb +   cc);
        varp += val8 * (   aa +   bb +   cc);
      }

      // normalize
      varp /= double( total_weight);

      // target resolution (2 sigma)
      const double sigma( m_Resolution / 2);

      linal::Vector< double> unknown_fac( s_MaxKernel, double( 0));
      unknown_fac( e_Gaussian)          =  0.0;
      unknown_fac( e_Triangular)        =  1.0;
      unknown_fac( e_SemiEpanechnikov)  =  1.5;
      unknown_fac( e_Epanechnikov)      =  2.0;
      unknown_fac( e_HardSphere)        = 60.0;

      linal::Vector< double> radius_half( s_MaxKernel, double( 0));
      radius_half( e_Gaussian)         = sigma * sqrt( log( 2.0)) / sqrt(1.5);
      radius_half( e_Triangular)       = sigma / ( exp( ( 1.0 / unknown_fac( e_Triangular))       * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_Triangular))       / (5.0*(5.0+unknown_fac( e_Triangular)))));
      radius_half( e_SemiEpanechnikov) = sigma / ( exp( ( 1.0 / unknown_fac( e_SemiEpanechnikov)) * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_SemiEpanechnikov)) / (5.0*(5.0+unknown_fac( e_SemiEpanechnikov)))));
      radius_half( e_Epanechnikov)     = sigma / ( exp( ( 1.0 / unknown_fac( e_Epanechnikov))     * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_Epanechnikov))     / (5.0*(5.0+unknown_fac( e_Epanechnikov)))));
      radius_half( e_HardSphere)       = sigma / ( exp( ( 1.0 / unknown_fac( e_HardSphere))       * log( 2.0)) * sqrt( 3.0*( 3.0+unknown_fac( e_HardSphere))       / (5.0*(5.0+unknown_fac( e_HardSphere)))));

      linal::Vector< double> radius_cut( s_MaxKernel, double( 0));
      radius_cut( e_Gaussian)         = sqrt( 3.0) * sigma;
      radius_cut( e_Triangular)       = ( exp( ( 1.0 / unknown_fac( e_Triangular))       * log( 2.0))) * radius_half( e_Triangular);
      radius_cut( e_SemiEpanechnikov) = ( exp( ( 1.0 / unknown_fac( e_SemiEpanechnikov)) * log( 2.0))) * radius_half( e_SemiEpanechnikov);
      radius_cut( e_Epanechnikov)     = ( exp( ( 1.0 / unknown_fac( e_Epanechnikov))     * log( 2.0))) * radius_half( e_Epanechnikov);
      radius_cut( e_HardSphere)       = ( exp( ( 1.0 / unknown_fac( e_HardSphere))       * log( 2.0))) * radius_half( e_HardSphere);

      Kernel kernel( m_SmoothingKernel);

      // if the half radius is too small than apply no smoothing
      if( radius_half( e_Gaussian) / m_GridSpacing.X() < 1.0)
      {
        kernel = e_NoSmoothing;
      }

      double kmsd( 0);
      if( kernel == e_Gaussian)
      {
        kmsd = math::Sqr( sigma) / ( m_GridSpacing.X() * m_GridSpacing.X());
      }
      else
      {
        kmsd = math::Sqr( sigma) / math::Sqr( m_GridSpacing.X());
      }

      double varmap( 0);
      if( corrmode == 1)
      {
        varmap -= varp; // variances are additive for uncorrelated samples
      }
      else // ( corrmode == 2)
      {
        varmap = kmsd;
      }

      if( varmap < 0)
      {
        // lattice smoothing exceeds kernel size
        return Map();
      }

      // maximal spanning of density coming from one atom
      int exth( 0);
      double sigmamap( 0);

      // compute lattice noise corrected kernel maps
      switch( kernel)
      {
        case e_NoSmoothing:
          break;

        case e_Gaussian:
        {
          sigmamap = sqrt( varmap / 3.0);   // sigma-1D
          exth = int( ceil( 3 * sigmamap)); // truncate at 3 sigma-1D == sqrt(3) sigma-3D
          break;
        }

        case e_Triangular:
        case e_SemiEpanechnikov:
        case e_Epanechnikov:
        case e_HardSphere:
        {
          exth = int( ceil( radius_cut( kernel) / m_GridSpacing.X()));
          break;
        }

        default:
          break;
      }

      // allocate kernel map
      math::Tensor< double> kernel_map( 2 * exth + 1, 2 * exth + 1, 2 * exth + 1, double( 0));

      // smooth map depending on kernel
      switch( kernel)
      {
        case e_NoSmoothing:
          break;

        case e_Gaussian:
        {
          // write Gaussian within 3 sigma-1D to map
          const double bvalue( -1.0 / ( 2.0 * math::Sqr( sigmamap)));
          const double cvalue( 9 * math::Sqr( sigmamap));
          for( size_t indz( 0); indz < kernel_map.NumberLayers(); ++indz)
          {
            for( size_t indy( 0); indy < kernel_map.GetNumberRows(); ++indy)
            {
              for( size_t indx( 0); indx < kernel_map.GetNumberCols(); ++indx)
              {
                const double dsqu( math::Sqr( indx - exth) + math::Sqr( indy - exth) + math::Sqr( indz - exth));
                if( dsqu < cvalue)
                {
                  kernel_map( indz, indy, indx) = kernel_amplitude * exp( dsqu * bvalue);
                }
              }
            }
          }

          break;
        }

        case e_Triangular:
        case e_SemiEpanechnikov:
        case e_Epanechnikov:
        case e_HardSphere:
        {
          // write kernel to map
          const double bvalue( 0.5 * exp( -unknown_fac( kernel) * log( radius_half( kernel))) * exp( unknown_fac( kernel) * log( m_GridSpacing.X())));
          for( size_t indz( 0); indz < kernel_map.NumberLayers(); ++indz)
          {
            for( size_t indy( 0); indy < kernel_map.GetNumberRows(); ++indy)
            {
              for( size_t indx( 0); indx < kernel_map.GetNumberCols(); ++indx)
              {
                const double dsqu
                (
                  exp
                  (
                    ( unknown_fac( kernel) / 2.0)
                    * std::log
                    (
                      double
                      (
                          math::Sqr( indx - exth)
                        + math::Sqr( indy - exth)
                        + math::Sqr( indz - exth)
                      )
                    )
                  )
                );
                kernel_map( indz, indy, indx) = std::max( double( 0.0), kernel_amplitude * ( 1.0 - dsqu * bvalue));
              }
            }
          }
          break;
        }

        default:
          break;
      }

      // convolve and write output
      switch( kernel)
      {
        case e_NoSmoothing:
        {
          // index of map
          const linal::VectorND< int, 3> index
          (
            int( floor( grid2.X() / m_GridSpacing.X() + 0.5)),
            int( floor( grid2.Y() / m_GridSpacing.Y() + 0.5)),
            int( floor( grid2.Z() / m_GridSpacing.Z() + 0.5))
          );

          // intervals
          linal::VectorND< int, 3> intervals;
          for( size_t i( 0); i < 3; ++i)
          {
            if( index( i) <= -int( ext_atom_map( i)))
            {
              intervals( i) = math::Absolute( index( i));
            }
            else if( index( i) >= 0)
            {
              intervals( i) = index( i) + ext_atom_map( i) - 1;
            }
            else
            {
              intervals( i) = ext_atom_map( i) - 1;
            }
          }

          // length
          const linal::Vector3D length
          (
            m_GridSpacing.X() * double( intervals( 0)),
            m_GridSpacing.Y() * double( intervals( 1)),
            m_GridSpacing.Z() * double( intervals( 2))
          );

//          BCL_Message
//          (
//            util::Message::e_Standard,
//            util::Format()( "kernel width should be large relative to lattice voxel spacing!\n") +
//            "Spatial resolution (2 sigma) from lattice smoothing alone: " +
//            util::Format().W( 6).FFP( 3)( m_GridSpacing.X() * sqrt(varp) * 2.0) + "\n" +
//            "If kernel smoothing desired, increase kernel width or decrease lattice spacing."
//          );

          return Map
          (
            atom_map,
            index,
            intervals,
            length,
            m_GridSpacing,
            Map::GetDefaultAngle(),
            Map::GetDefaultAxis(),
            linal::Vector3D( 0.0)
          );
        }
        case e_Gaussian:
        case e_Triangular:
        case e_SemiEpanechnikov:
        case e_Epanechnikov:
        case e_HardSphere:
        {
          // Convolving lattice with kernel

          // allocate output density map
          const storage::VectorND< 3, size_t> ext_output // extension output map
          (
            size_t( ceil( ( maxcoord.X() - mincoord.X()) / m_GridSpacing.X()) + 2 * exth + 2 * m_Margin + 1),
            size_t( ceil( ( maxcoord.Y() - mincoord.Y()) / m_GridSpacing.Y()) + 2 * exth + 2 * m_Margin + 1),
            size_t( ceil( ( maxcoord.Z() - mincoord.Z()) / m_GridSpacing.Z()) + 2 * exth + 2 * m_Margin + 1)
          );
          const linal::Vector3D grid3
          (
            mincoord.X() - ( exth + m_Margin) * m_GridSpacing.X(),
            mincoord.Y() - ( exth + m_Margin) * m_GridSpacing.Y(),
            mincoord.Z() - ( exth + m_Margin) * m_GridSpacing.Z()
          );
          math::Tensor< double> output_map( ext_output.Third(), ext_output.Second(), ext_output.First(), double( 0));
          const size_t extxy2( ext_atom_map.First() * ext_atom_map.Second());
          size_t count( 0);
          for
          (
            const double *density_val( atom_map.Begin()), *density_val_end( atom_map.End());
            density_val != density_val_end;
            ++density_val, ++count
          )
          {
            if( *density_val != 0.0)
            {
              int indv = int( count);
              int k = indv / extxy2;
              indv -= k * extxy2;
              int j = indv / ext_atom_map.First();
              int i = indv - j * ext_atom_map.First();
              for( size_t indz( 0); indz < kernel_map.NumberLayers(); ++indz)
              {
                for( size_t indy( 0); indy < kernel_map.GetNumberRows(); ++indy)
                {
                  for( size_t indx( 0); indx < kernel_map.GetNumberCols(); ++indx)
                  {
                    output_map( k + indz, j + indy, i + indx) += kernel_map( indz, indy, indx) * ( *density_val);
                  }
                }
              }
            }
          }

//          if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Verbose))
//          {
//            std::ostringstream msg;
//            msg << "Lattice smoothing (sigma = atom rmsd) done: X Y Z "
//                << m_GridSpacing.X() * sqrt(varp) << ' ' << m_GridSpacing.Y() * sqrt(varp) << ' ' << m_GridSpacing.Z() * sqrt(varp)
//                << " Angstroem" << '\n'
//                << "Spatial resolution (2 sigma) of output map: ";
//            if( corrmode == 2)
//            {
//              const double reso( 2 * sqrt( math::Sqr( sigma) + varp * m_GridSpacing.X() * m_GridSpacing.X())); // variances are additive for uncorrelated samples
//              msg << reso << " slightly larger than target resolution due to uncorrected lattice smoothing";
//            }
//            else
//            {
//              msg << m_Resolution;
//            }
//            BCL_MessageVrb( msg.str());
//          }

          // index of map
          const linal::VectorND< int, 3> index
          (
            int( floor( grid3.X() / m_GridSpacing.X() + 0.5)),
            int( floor( grid3.Y() / m_GridSpacing.Y() + 0.5)),
            int( floor( grid3.Z() / m_GridSpacing.Z() + 0.5))
          );

          // intervals
          linal::VectorND< int, 3> intervals;
          for( size_t i( 0); i < 3; ++i)
          {
            if( index( i) <= -int( ext_output( i)))
            {
              intervals( i) = math::Absolute( index( i));
            }
            else if( index( i) >= 0)
            {
              intervals( i) = index( i) + ext_output( i) - 1;
            }
            else
            {
              intervals( i) = ext_output( i) - 1;
            }
          }

          // length
          const linal::Vector3D length
          (
            m_GridSpacing.X() * double( intervals( 0)),
            m_GridSpacing.Y() * double( intervals( 1)),
            m_GridSpacing.Z() * double( intervals( 2))
          );

          // construct map and return
          return Map
            (
              output_map,
              index,
              intervals,
              length,
              m_GridSpacing,
              Map::GetDefaultAngle(),
              Map::GetDefaultAxis(),
              linal::Vector3D( 0.0)
            );
        }

        default:
        {
          return Map();
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SimulateDefault::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_GridSpacing    , ISTREAM);
      io::Serialize::Read( m_Resolution     , ISTREAM);
      io::Serialize::Read( m_Margin         , ISTREAM);
      io::Serialize::Read( m_SmoothingKernel, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SimulateDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_GridSpacing    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Resolution     , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Margin         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SmoothingKernel, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace density
} // namespace bcl
