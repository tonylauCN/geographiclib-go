package geodesic

import (
	"math"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"
)

const NumIt int = 10

var Eps = 0.01 * math.Sqrt(epsilon)

type GnomonicData struct {
	Lat0 float64 // latitude of center point of projection (degrees).
	Lon0 float64 // longitude of center point of projection (degrees).
	Lat  float64 // latitude of point (degrees).
	Lon  float64 // longitude of point (degrees).
	X    float64 // easting of point (meters).
	Y    float64 // northing of point (meters).
	Azi  float64 // azimuth of geodesic at point (degrees).
	Rk   float64 // reciprocal of azimuthal scale at point.
}

type Gnomonic struct {
	earth *Geodesic
	a     float64
	f     float64
}

func NewGnomonic(geo *Geodesic) *Gnomonic {
	g := &Gnomonic{
		earth: geo,
		a:     geo.EquatorialRadius(),
		f:     geo.Flattening(),
	}
	return g
}

func EarthGnomonic() *Gnomonic {

	return &Gnomonic{
		earth: WGS84,
		a:     WGS84_a,
		f:     WGS84_f,
	}
}

// Forward projection, from geographic to gnomonic.
func (g *Gnomonic) Forward(lat0, lon0, lat, lon float64) GnomonicData {
	inv := g.earth.InverseWithCapabilities(lat0, lon0, lat, lon,
		capabilities.Azimuth|capabilities.GeodesicScale|capabilities.ReducedLength)
	fwd := GnomonicData{lat0, lon0, lat, lon, math.NaN(), math.NaN(), inv.Azi2, inv.M12}
	if inv.M12 > 0.0 {
		rho := inv.M12Reduced / inv.M12
		first, second := sincosd(inv.Azi1)
		fwd.X = rho * first
		fwd.Y = rho * second
	}

	return fwd
}

// Reverse projection, from gnomonic to geographic.
func (g *Gnomonic) Reverse(lat0, lon0, x, y float64) GnomonicData {
	rev := GnomonicData{lat0, lon0, math.NaN(), math.NaN(), x, y, math.NaN(), math.NaN()}
	azi0 := atan2d(x, y)
	rho := math.Hypot(x, y)
	s := g.a * math.Atan(rho/g.a)
	little := rho <= g.a
	if !little {
		rho = 1.0 / rho
	}

	line := g.earth.LineWithCapabilities(lat0, lon0, azi0, capabilities.Latitude|capabilities.Longitude|capabilities.Azimuth|
		capabilities.DistanceIn|capabilities.ReducedLength|capabilities.GeodesicScale)
	count := NumIt
	trip := 0
	var pos Data

	for ; count > 0; count-- {
		pos = line.PositionWithCapabilities(s, capabilities.Latitude|capabilities.Longitude|capabilities.Azimuth|
			capabilities.DistanceIn|capabilities.ReducedLength|capabilities.GeodesicScale)
		if trip > 0 {
			break
		}

		ds := 0.0
		if little {
			ds = (pos.M12Reduced/pos.M12 - rho) * pos.M12 * pos.M12
		} else {
			ds = (rho - pos.M12/pos.M12Reduced) * pos.M12Reduced * pos.M12Reduced
		}
		s -= ds
		if math.Abs(ds) <= Eps*g.a {
			trip++
		}
	}

	if trip == 0 {
		return rev
	} else {
		rev.Lat = pos.Lat2
		rev.Lon = pos.Lon2
		rev.Azi = pos.Azi2
		rev.Rk = pos.M12
		return rev
	}
}
