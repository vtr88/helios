#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

/* Altere somente este bloco para usar outra cidade. */
#define DEFAULT_LOCATION_NAME "Tanabi, SP"
#define DEFAULT_LOCATION_LATITUDE (-20.63)
#define DEFAULT_LOCATION_LONGITUDE (-49.65)
#define DEFAULT_LOCATION_TIMEZONE "America/Sao_Paulo"

#define PI 3.14159265358979323846
#define RAD (PI / 180.0)
#define DAY_SECONDS 86400.0
#define J1970 2440588.0
#define J2000 2451545.0
#define J0 0.0009
#define EARTH_RADIUS_KM 6378.14

static const char *location_name = DEFAULT_LOCATION_NAME;
static double location_latitude = DEFAULT_LOCATION_LATITUDE;
static double location_longitude = DEFAULT_LOCATION_LONGITUDE;
static char location_timezone[128] = DEFAULT_LOCATION_TIMEZONE;

#define COLOR_BOLD "\033[1m"
#define COLOR_ORANGE "\033[38;5;208m"
#define COLOR_OLIVE "\033[38;5;142m"
#define COLOR_GOLD "\033[38;5;220m"
#define COLOR_DIM "\033[38;5;245m"
#define COLOR_RESET "\033[0m"

static bool color_output;

static const char *color(const char *sequence)
{
	return color_output ? sequence : "";
}

typedef struct { double ra, dec, dist; } Coordinates;
typedef struct { double altitude, azimuth, distance; } Position;
typedef struct { double fraction, phase; bool waxing; } Illumination;

enum SunEvent {
	NADIR, NIGHT_END, NAUTICAL_DAWN, DAWN, SUNRISE, SUNRISE_END,
	GOLDEN_HOUR_END, SOLAR_NOON, GOLDEN_HOUR, SUNSET_START, SUNSET,
	DUSK, NAUTICAL_DUSK, NIGHT, SUN_EVENT_COUNT
};

typedef struct { double at[SUN_EVENT_COUNT]; } SunTimes;
typedef struct { double rise, set; bool always_up, always_down; } MoonTimes;

static double to_days(double unix_seconds)
{
	return unix_seconds / DAY_SECONDS - 0.5 + J1970 - J2000;
}

static double from_julian(double julian)
{
	return (julian + 0.5 - J1970) * DAY_SECONDS;
}

static double delta_t(double d)
{
	double y = 2000.0 + d / 365.2425;
	double t;
	if (y < 1920.0) {
		t = y - 1900.0;
		return -2.79 + t * (1.494119 + t * (-0.0598939 + t * (0.0061966 - t * 0.000197)));
	}
	if (y < 1941.0) {
		t = y - 1920.0;
		return 21.20 + t * (0.84493 + t * (-0.076100 + t * 0.0020936));
	}
	if (y < 1961.0) {
		t = y - 1950.0;
		return 29.07 + t * (0.407 + t * (-1.0 / 233.0 + t / 2547.0));
	}
	if (y < 1986.0) {
		t = y - 1975.0;
		return 45.45 + t * (1.067 + t * (-1.0 / 260.0 - t / 718.0));
	}
	if (y < 2005.0) {
		t = y - 2000.0;
		return 63.86 + t * (0.3345 + t * (-0.060374 + t * (0.0017275 + t * (0.000651814 + t * 0.00002373599))));
	}
	if (y < 2050.0) {
		t = y - 2000.0;
		return 62.92 + t * (0.32217 + t * 0.005589);
	}
	t = (y - 1820.0) / 100.0;
	return -20.0 + 32.0 * t * t - 0.5628 * (2150.0 - y);
}

static double to_days_tt(double d)
{
	return d + delta_t(d) / DAY_SECONDS;
}

static double right_ascension(double l, double b)
{
	const double e = RAD * 23.4397;
	return atan2(sin(l) * cos(e) - tan(b) * sin(e), cos(l));
}

static double declination(double l, double b)
{
	const double e = RAD * 23.4397;
	return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l));
}

static double altitude(double hour_angle, double latitude, double dec)
{
	return asin(sin(latitude) * sin(dec) + cos(latitude) * cos(dec) * cos(hour_angle));
}

static double azimuth(double hour_angle, double latitude, double dec)
{
	double degrees = atan2(sin(hour_angle), cos(hour_angle) * sin(latitude) - tan(dec) * cos(latitude)) / RAD;
	degrees = fmod(degrees + 540.0, 360.0);
	return degrees < 0.0 ? degrees + 360.0 : degrees;
}

static double sidereal_time(double d, double lw)
{
	return RAD * (280.46061837 + 360.98564736629 * d) - lw;
}

static double astro_refraction(double h)
{
	if (h < 0.0) h = 0.0;
	return 0.0002967 / tan(h + 0.00312536 / (h + 0.08901179));
}

static Coordinates sun_coords(double d)
{
	double t = d / 36525.0;
	double l0 = RAD * (280.46646 + t * (36000.76983 + t * 0.0003032));
	double m = RAD * (357.52911 + t * (35999.05029 - t * 0.0001537));
	double sin_m = sin(m), cos_m = cos(m);
	double c = RAD * ((1.914602 - t * (0.004817 + t * 0.000014)) * sin_m
		+ (0.019993 - 0.000101 * t) * 2.0 * sin_m * cos_m
		+ 0.000289 * sin_m * (3.0 - 4.0 * sin_m * sin_m));
	double om = RAD * (125.04 - 1934.136 * t);
	double l = l0 + c - RAD * (0.00569 + 0.00478 * sin(om));
	double e = RAD * (23.439291 - t * (0.0130042 + t * (0.00000016 - t * 0.000000504)))
		+ RAD * 0.00256 * cos(om);
	Coordinates result = {
		.ra = atan2(cos(e) * sin(l), cos(l)),
		.dec = asin(sin(e) * sin(l)),
		.dist = 149598000.0
	};
	return result;
}

static Position sun_position(double unix_seconds)
{
	double lw = RAD * -location_longitude;
	double phi = RAD * location_latitude;
	double d = to_days(unix_seconds);
	Coordinates c = sun_coords(to_days_tt(d));
	double hour_angle = sidereal_time(d, lw) - c.ra;
	double h = altitude(hour_angle, phi, c.dec);
	Position result = {
		.altitude = (h + astro_refraction(h)) / RAD,
		.azimuth = azimuth(hour_angle, phi, c.dec),
		.distance = c.dist
	};
	return result;
}

static double wrap_pi(double angle)
{
	return angle - 2.0 * PI * round(angle / (2.0 * PI));
}

static double solar_transit(double dt, double lw)
{
	for (int i = 0; i < 3; i++) {
		double h = wrap_pi(sidereal_time(dt, lw) - sun_coords(to_days_tt(dt)).ra);
		dt -= h / (2.0 * PI);
	}
	return dt;
}

static double sun_set_j(double h0, double dt, int sign, double lw, double phi, double dec_at_transit)
{
	double cos_h0 = (sin(h0) - sin(phi) * sin(dec_at_transit)) / (cos(phi) * cos(dec_at_transit));
	if (cos_h0 < -1.0 || cos_h0 > 1.0) return NAN;
	double d = dt + sign * acos(cos_h0) / (2.0 * PI);
	for (int i = 0; i < 2; i++) {
		Coordinates c = sun_coords(to_days_tt(d));
		double hour_angle = wrap_pi(sidereal_time(d, lw) - c.ra);
		double h = altitude(hour_angle, phi, c.dec);
		double sin_h = cos(phi) * cos(c.dec) * sin(hour_angle);
		if (fabs(sin_h) < 1e-6) break;
		d += (h - h0) / (2.0 * PI * sin_h);
	}
	return d;
}

static SunTimes sun_times(double unix_seconds)
{
	static const struct { double angle; enum SunEvent morning, evening; } crossings[] = {
		{-0.833, SUNRISE, SUNSET}, {-0.3, SUNRISE_END, SUNSET_START},
		{-6.0, DAWN, DUSK}, {-12.0, NAUTICAL_DAWN, NAUTICAL_DUSK},
		{-18.0, NIGHT_END, NIGHT}, {6.0, GOLDEN_HOUR_END, GOLDEN_HOUR}
	};
	SunTimes result;
	for (int i = 0; i < SUN_EVENT_COUNT; i++) result.at[i] = NAN;
	double lw = RAD * -location_longitude;
	double phi = RAD * location_latitude;
	double d = round(round(to_days(unix_seconds)) - J0 - lw / (2.0 * PI));
	double transit = solar_transit(d + J0 + lw / (2.0 * PI), lw);
	double dec = sun_coords(to_days_tt(transit)).dec;
	result.at[SOLAR_NOON] = from_julian(transit + J2000);
	result.at[NADIR] = from_julian(transit + J2000 - 0.5);
	for (size_t i = 0; i < sizeof(crossings) / sizeof(crossings[0]); i++) {
		double rise = sun_set_j(crossings[i].angle * RAD, transit, -1, lw, phi, dec);
		double set = sun_set_j(crossings[i].angle * RAD, transit, 1, lw, phi, dec);
		if (!isnan(rise)) result.at[crossings[i].morning] = from_julian(rise + J2000);
		if (!isnan(set)) result.at[crossings[i].evening] = from_julian(set + J2000);
	}
	return result;
}

/* Série lunar compacta usada pelo SunCalc clássico (Meeus, cap. 47). */
static Coordinates moon_coords(double d)
{
	double l = RAD * (218.316 + 13.176396 * d);
	double m = RAD * (134.963 + 13.064993 * d);
	double f = RAD * (93.272 + 13.229350 * d);
	double longitude = l + RAD * 6.289 * sin(m);
	double latitude = RAD * 5.128 * sin(f);
	Coordinates result = {
		.ra = right_ascension(longitude, latitude),
		.dec = declination(longitude, latitude),
		.dist = 385001.0 - 20905.0 * cos(m)
	};
	return result;
}

static Position moon_position(double unix_seconds)
{
	double lw = RAD * -location_longitude;
	double phi = RAD * location_latitude;
	double d = to_days(unix_seconds);
	Coordinates c = moon_coords(to_days_tt(d));
	double hour_angle = sidereal_time(d, lw) - c.ra;
	double geocentric_h = altitude(hour_angle, phi, c.dec);
	double h = geocentric_h - asin(EARTH_RADIUS_KM / c.dist * cos(geocentric_h));
	Position result = {
		.altitude = (h + astro_refraction(h)) / RAD,
		.azimuth = azimuth(hour_angle, phi, c.dec),
		.distance = c.dist
	};
	return result;
}

static Illumination moon_illumination(double unix_seconds)
{
	double d = to_days_tt(to_days(unix_seconds));
	Coordinates s = sun_coords(d), m = moon_coords(d);
	double phi = acos(sin(s.dec) * sin(m.dec) + cos(s.dec) * cos(m.dec) * cos(s.ra - m.ra));
	double inc = atan2(s.dist * sin(phi), m.dist - s.dist * cos(phi));
	double angle = atan2(cos(s.dec) * sin(s.ra - m.ra),
		sin(s.dec) * cos(m.dec) - cos(s.dec) * sin(m.dec) * cos(s.ra - m.ra));
	bool waxing = angle < 0.0;
	Illumination result = {
		.fraction = (1.0 + cos(inc)) / 2.0,
		.phase = 0.5 + 0.5 * inc * (waxing ? -1.0 : 1.0) / PI,
		.waxing = waxing
	};
	return result;
}

static double moon_height(double unix_seconds)
{
	Position p = moon_position(unix_seconds);
	return p.altitude + 0.2725 * asin(EARTH_RADIUS_KM / p.distance) / RAD + 0.09;
}

static double refine_moon_crossing(double unix_seconds)
{
	for (int i = 0; i < 2; i++) {
		double h = moon_height(unix_seconds);
		double dh = (moon_height(unix_seconds + 30.0) - moon_height(unix_seconds - 30.0)) / 60.0;
		if (fabs(dh) < 1e-9) break;
		unix_seconds -= h / dh;
	}
	return unix_seconds;
}

static MoonTimes moon_times(double local_midnight)
{
	MoonTimes result = {.rise = NAN, .set = NAN, .always_up = false, .always_down = false};
	double h0 = moon_height(local_midnight), h_max = h0;
	for (int i = 1; i <= 24; i += 2) {
		double h1 = moon_height(local_midnight + i * 3600.0);
		double h2 = moon_height(local_midnight + (i + 1) * 3600.0);
		if (h1 > h_max) h_max = h1;
		if (h2 > h_max) h_max = h2;
		double a = (h0 + h2) / 2.0 - h1;
		double b = (h2 - h0) / 2.0;
		if (fabs(a) < 1e-12) { h0 = h2; continue; }
		double xe = -b / (2.0 * a);
		double discriminant = b * b - 4.0 * a * h1;
		double ye = (a * xe + b) * xe + h1;
		int roots = 0;
		double x1 = 0.0, x2 = 0.0;
		if (discriminant >= 0.0) {
			double dx = sqrt(discriminant) / (fabs(a) * 2.0);
			x1 = xe - dx;
			x2 = xe + dx;
			if (fabs(x1) <= 1.0) roots++;
			if (fabs(x2) <= 1.0) roots++;
			if (x1 < -1.0) x1 = x2;
		}
		if (roots == 1) {
			double crossing = local_midnight + (i + x1) * 3600.0;
			if (h0 < 0.0) result.rise = crossing; else result.set = crossing;
		} else if (roots == 2) {
			result.rise = local_midnight + (i + (ye < 0.0 ? x2 : x1)) * 3600.0;
			result.set = local_midnight + (i + (ye < 0.0 ? x1 : x2)) * 3600.0;
		}
		if (!isnan(result.rise) && !isnan(result.set)) break;
		h0 = h2;
	}
	if (!isnan(result.rise)) result.rise = refine_moon_crossing(result.rise);
	if (!isnan(result.set)) result.set = refine_moon_crossing(result.set);
	if (isnan(result.rise) && isnan(result.set)) {
		result.always_up = h_max > 0.0;
		result.always_down = h_max <= 0.0;
	}
	return result;
}

static const char *direction(double degrees)
{
	static const char *names[] = {
		"N", "NNE", "NE", "ENE", "L", "ESE", "SE", "SSE",
		"S", "SSO", "SO", "OSO", "O", "ONO", "NO", "NNO"
	};
	return names[(int)floor((degrees + 11.25) / 22.5) % 16];
}

static const char *phase_name(double phase)
{
	static const char *names[] = {
		"nova", "crescente", "quarto crescente", "gibosa crescente",
		"cheia", "gibosa minguante", "quarto minguante", "minguante"
	};
	return names[(int)floor(phase * 8.0 + 0.5) % 8];
}

static void print_time_value(double unix_seconds)
{
	if (isnan(unix_seconds)) { printf("--:--"); return; }
	time_t value = (time_t)llround(unix_seconds);
	struct tm local;
	char buffer[16];
	localtime_r(&value, &local);
	strftime(buffer, sizeof(buffer), "%H:%M", &local);
	printf("%s", buffer);
}

static void print_event(const char *name, double unix_seconds, const char *explanation, const char *label_color)
{
	printf("  %s", color(COLOR_GOLD));
	print_time_value(unix_seconds);
	printf("%s  %s%s%s", color(COLOR_RESET), color(label_color), name, color(COLOR_RESET));
	if (explanation) printf(" %s— %s%s", color(COLOR_DIM), explanation, color(COLOR_RESET));
	putchar('\n');
}

static bool parse_date(const char *text, struct tm *date)
{
	int year, month, day;
	char tail;
	if (strlen(text) != 10 || text[4] != '-' || text[7] != '-') return false;
	if (sscanf(text, "%d-%d-%d%c", &year, &month, &day, &tail) != 3) return false;
	if (year < 1900 || year > 2150 || month < 1 || month > 12 || day < 1 || day > 31) return false;
	memset(date, 0, sizeof(*date));
	date->tm_year = year - 1900;
	date->tm_mon = month - 1;
	date->tm_mday = day;
	date->tm_isdst = -1;
	time_t checked = mktime(date);
	return checked != (time_t)-1 && date->tm_year == year - 1900
		&& date->tm_mon == month - 1 && date->tm_mday == day;
}

static bool timezone_exists(const char *timezone)
{
	if (!timezone[0] || timezone[0] == '/' || strstr(timezone, "..")) return false;
	for (const unsigned char *p = (const unsigned char *)timezone; *p; p++) {
		if (!(('A' <= *p && *p <= 'Z') || ('a' <= *p && *p <= 'z')
				|| ('0' <= *p && *p <= '9') || *p == '/' || *p == '_'
				|| *p == '-' || *p == '+')) return false;
	}
	char path[256];
	if (snprintf(path, sizeof(path), "/usr/share/zoneinfo/%s", timezone) >= (int)sizeof(path)) return false;
	return access(path, R_OK) == 0;
}

static bool parse_local(const char *spec)
{
	char buffer[512];
	if (strlen(spec) >= sizeof(buffer)) return false;
	strcpy(buffer, spec);
	char *first = strchr(buffer, ',');
	if (!first) return false;
	*first++ = '\0';
	char *second = strchr(first, ',');
	if (!second || strchr(second + 1, ',')) return false;
	*second++ = '\0';

	char *end;
	errno = 0;
	double latitude = strtod(buffer, &end);
	if (errno || *end || latitude < -90.0 || latitude > 90.0) return false;
	errno = 0;
	double longitude = strtod(first, &end);
	if (errno || *end || longitude < -180.0 || longitude > 180.0) return false;
	if (!timezone_exists(second) || strlen(second) >= sizeof(location_timezone)) return false;

	location_name = "Local temporário";
	location_latitude = latitude;
	location_longitude = longitude;
	strcpy(location_timezone, second);
	return true;
}

static int hex_digit(unsigned char c)
{
	if ('0' <= c && c <= '9') return c - '0';
	if ('a' <= c && c <= 'f') return c - 'a' + 10;
	if ('A' <= c && c <= 'F') return c - 'A' + 10;
	return -1;
}

static bool append_utf8(char *out, size_t capacity, size_t *length, unsigned codepoint)
{
	unsigned char bytes[3];
	size_t count;
	if (codepoint <= 0x7f) { bytes[0] = codepoint; count = 1; }
	else if (codepoint <= 0x7ff) {
		bytes[0] = 0xc0 | (codepoint >> 6);
		bytes[1] = 0x80 | (codepoint & 0x3f);
		count = 2;
	} else {
		bytes[0] = 0xe0 | (codepoint >> 12);
		bytes[1] = 0x80 | ((codepoint >> 6) & 0x3f);
		bytes[2] = 0x80 | (codepoint & 0x3f);
		count = 3;
	}
	if (*length + count >= capacity) return false;
	for (size_t i = 0; i < count; i++) out[(*length)++] = bytes[i];
	return true;
}

static const char *json_field(const char *object, const char *key)
{
	char pattern[80];
	if (snprintf(pattern, sizeof(pattern), "\"%s\"", key) >= (int)sizeof(pattern)) return NULL;
	const char *p = strstr(object, pattern);
	if (!p) return NULL;
	p += strlen(pattern);
	while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
	if (*p++ != ':') return NULL;
	while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
	return p;
}

static bool json_string(const char *object, const char *key, char *out, size_t capacity)
{
	const unsigned char *p = (const unsigned char *)json_field(object, key);
	if (!p || *p++ != '"') return false;
	size_t length = 0;
	while (*p && *p != '"') {
		unsigned value = *p++;
		if (value == '\\') {
			value = *p++;
			if (value == 'n') value = '\n';
			else if (value == 'r') value = '\r';
			else if (value == 't') value = '\t';
			else if (value == 'b') value = '\b';
			else if (value == 'f') value = '\f';
			else if (value == 'u') {
				unsigned codepoint = 0;
				for (int i = 0; i < 4; i++) {
					int digit = hex_digit(*p++);
					if (digit < 0) return false;
					codepoint = codepoint * 16 + (unsigned)digit;
				}
				if (!append_utf8(out, capacity, &length, codepoint)) return false;
				continue;
			}
		}
		if (length + 1 >= capacity) return false;
		out[length++] = (char)value;
	}
	if (*p != '"') return false;
	out[length] = '\0';
	return true;
}

static bool json_number(const char *object, const char *key, double *value)
{
	const char *p = json_field(object, key);
	if (!p) return false;
	char *end;
	errno = 0;
	*value = strtod(p, &end);
	return !errno && end != p;
}

static char *fetch_locations(const char *query)
{
	int descriptors[2];
	if (pipe(descriptors) != 0) return NULL;
	pid_t child = fork();
	if (child < 0) {
		close(descriptors[0]); close(descriptors[1]);
		return NULL;
	}
	if (child == 0) {
		close(descriptors[0]);
		if (dup2(descriptors[1], STDOUT_FILENO) < 0) _exit(127);
		close(descriptors[1]);
		char *name = malloc(strlen(query) + 6);
		if (!name) _exit(127);
		sprintf(name, "name=%s", query);
		execlp("curl", "curl", "-fsS", "--max-time", "10", "--get",
			"--data-urlencode", name, "--data", "count=5", "--data", "language=pt",
			"--data", "format=json", "--user-agent", "helios/0.1",
			"https://geocoding-api.open-meteo.com/v1/search", (char *)NULL);
		_exit(127);
	}

	close(descriptors[1]);
	size_t length = 0, capacity = 4096;
	char *response = malloc(capacity);
	if (!response) { close(descriptors[0]); waitpid(child, NULL, 0); return NULL; }
	for (;;) {
		if (length + 2048 + 1 > capacity) {
			if (capacity >= 1024 * 1024) { free(response); response = NULL; break; }
			capacity *= 2;
			char *larger = realloc(response, capacity);
			if (!larger) { free(response); response = NULL; break; }
			response = larger;
		}
		ssize_t count = read(descriptors[0], response + length, capacity - length - 1);
		if (count == 0) break;
		if (count < 0) { free(response); response = NULL; break; }
		length += (size_t)count;
	}
	close(descriptors[0]);
	int status;
	if (waitpid(child, &status, 0) < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0) {
		free(response);
		return NULL;
	}
	if (response) response[length] = '\0';
	return response;
}

static int find_location(const char *query, const char *program)
{
	if (strlen(query) < 2) {
		fprintf(stderr, "busca muito curta: digite ao menos dois caracteres\n");
		return 2;
	}
	char *json = fetch_locations(query);
	if (!json) {
		fprintf(stderr, "não foi possível consultar locais; verifique a rede e se curl está instalado\n");
		return 1;
	}
	const char *results = strstr(json, "\"results\"");
	if (!results || !(results = strchr(results, '['))) {
		printf("nenhum local encontrado para \"%s\"\n", query);
		free(json);
		return 1;
	}
	printf("Resultados para \"%s\":\n", query);
	int found = 0;
	const char *p = results + 1;
	while (*p && *p != ']' && found < 5) {
		while (*p && *p != '{' && *p != ']') p++;
		if (*p != '{') break;
		const char *start = p;
		int depth = 0;
		bool in_string = false, escaped = false;
		do {
			if (in_string) {
				if (escaped) escaped = false;
				else if (*p == '\\') escaped = true;
				else if (*p == '"') in_string = false;
			} else {
				if (*p == '"') in_string = true;
				else if (*p == '{') depth++;
				else if (*p == '}') depth--;
			}
			p++;
		} while (*p && depth > 0);
		if (depth != 0) break;
		size_t size = (size_t)(p - start);
		char *object = malloc(size + 1);
		if (!object) break;
		memcpy(object, start, size);
		object[size] = '\0';
		char name[160], admin[160] = "", country[160] = "", timezone[128];
		double latitude, longitude;
		if (json_string(object, "name", name, sizeof(name))
				&& json_string(object, "timezone", timezone, sizeof(timezone))
				&& json_number(object, "latitude", &latitude)
				&& json_number(object, "longitude", &longitude)) {
			json_string(object, "admin1", admin, sizeof(admin));
			json_string(object, "country", country, sizeof(country));
			found++;
			printf("\n%d. %s", found, name);
			if (admin[0] && strcmp(admin, name) != 0) printf(", %s", admin);
			if (country[0]) printf(", %s", country);
			printf("\n   latitude %.5f · longitude %.5f · fuso %s\n", latitude, longitude, timezone);
			printf("   %s --local=%.5f,%.5f,%s\n", program, latitude, longitude, timezone);
		}
		free(object);
	}
	free(json);
	if (!found) {
		printf("nenhum local encontrado para \"%s\"\n", query);
		return 1;
	}
	return 0;
}

static void usage(const char *program)
{
	printf("uso: %s [--local=LAT,LONG,FUSO] [AAAA-MM-DD]\n", program);
	printf("     %s --ache \"cidade[, estado ou país]\"\n", program);
	printf("sem data, usa o instante e o dia civil atuais do local escolhido.\n");
}

int main(int argc, char **argv)
{
	const char *term = getenv("TERM");
	color_output = isatty(STDOUT_FILENO) && getenv("NO_COLOR") == NULL
		&& (!term || strcmp(term, "dumb") != 0);
	if (argc >= 2 && (strcmp(argv[1], "--ache") == 0 || strncmp(argv[1], "--ache=", 8) == 0)) {
		const char *query;
		if (strncmp(argv[1], "--ache=", 8) == 0) {
			if (argc != 2 || !argv[1][8]) { usage(argv[0]); return 2; }
			query = argv[1] + 8;
		} else {
			if (argc != 3) { usage(argv[0]); return 2; }
			query = argv[2];
		}
		return find_location(query, argv[0]);
	}

	const char *date_arg = NULL;
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			usage(argv[0]);
			return 0;
		}
		const char *local_spec = NULL;
		if (strncmp(argv[i], "--local=", 8) == 0) local_spec = argv[i] + 8;
		else if (strcmp(argv[i], "--local") == 0 && i + 1 < argc) local_spec = argv[++i];
		if (local_spec) {
			if (!parse_local(local_spec)) {
				fprintf(stderr, "local inválido: use LATITUDE,LONGITUDE,FUSO_IANA\n");
				return 2;
			}
			continue;
		}
		if (argv[i][0] == '-' || date_arg) { usage(argv[0]); return 2; }
		date_arg = argv[i];
	}

	if (setenv("TZ", location_timezone, 1) != 0) {
		fprintf(stderr, "erro ao definir o fuso: %s\n", strerror(errno));
		return 1;
	}
	tzset();
	time_t now = time(NULL);
	struct tm date;
	bool is_today = date_arg == NULL;
	if (is_today) localtime_r(&now, &date);
	else if (!parse_date(date_arg, &date)) {
		fprintf(stderr, "data inválida: %s (use AAAA-MM-DD)\n", date_arg);
		return 2;
	}

	struct tm midnight_tm = date;
	midnight_tm.tm_hour = midnight_tm.tm_min = midnight_tm.tm_sec = 0;
	midnight_tm.tm_isdst = -1;
	time_t midnight = mktime(&midnight_tm);
	struct tm noon_tm = date;
	noon_tm.tm_hour = 12;
	noon_tm.tm_min = noon_tm.tm_sec = 0;
	noon_tm.tm_isdst = -1;
	time_t noon = mktime(&noon_tm);
	double observation = is_today ? (double)now : (double)noon;

	char date_buffer[32], time_buffer[32];
	struct tm observation_tm;
	time_t observation_time = (time_t)observation;
	localtime_r(&observation_time, &observation_tm);
	strftime(date_buffer, sizeof(date_buffer), "%d/%m/%Y", &observation_tm);
	strftime(time_buffer, sizeof(time_buffer), "%H:%M %Z", &observation_tm);

	Position sun = sun_position(observation), moon = moon_position(observation);
	Illumination light = moon_illumination(observation);
	SunTimes st = sun_times((double)noon);
	MoonTimes mt = moon_times((double)midnight);

	printf("%s%s%s — %s%s\n", color(COLOR_BOLD), color(COLOR_ORANGE), location_name, date_buffer, color(COLOR_RESET));
	printf("%s%.2f, %.2f · %s%s\n", color(COLOR_DIM), location_latitude, location_longitude, location_timezone, color(COLOR_RESET));
	printf("\n%s%sAGORA · %s%s\n", color(COLOR_BOLD), color(COLOR_GOLD), time_buffer, color(COLOR_RESET));
	printf("  %sSol%s  altura %5.1f° · direção %s (%5.1f°)\n",
		color(COLOR_ORANGE), color(COLOR_RESET), sun.altitude, direction(sun.azimuth), sun.azimuth);
	printf("  %sLua%s  altura %5.1f° · direção %s (%5.1f°) · %.0f km\n",
		color(COLOR_OLIVE), color(COLOR_RESET), moon.altitude, direction(moon.azimuth), moon.azimuth, moon.distance);
	printf("       %s · %.0f%% iluminada · %s\n",
		phase_name(light.phase), light.fraction * 100.0, light.waxing ? "crescendo" : "minguando");
	printf("  %saltura: 0° é o horizonte · direção: pontos da bússola%s\n", color(COLOR_DIM), color(COLOR_RESET));

	printf("\n%s%sSOL%s\n", color(COLOR_BOLD), color(COLOR_ORANGE), color(COLOR_RESET));
	print_event("Noite astronômica termina", st.at[NIGHT_END], "primeira claridade", COLOR_ORANGE);
	print_event("Aurora náutica", st.at[NAUTICAL_DAWN], "o horizonte começa a aparecer", COLOR_ORANGE);
	print_event("Aurora civil", st.at[DAWN], "já há luz para atividades externas", COLOR_ORANGE);
	print_event("Nascer do sol", st.at[SUNRISE], "a primeira borda aparece", COLOR_ORANGE);
	print_event("Sol totalmente visível", st.at[SUNRISE_END], NULL, COLOR_ORANGE);
	print_event("Fim da hora dourada", st.at[GOLDEN_HOUR_END], "termina a luz suave da manhã", COLOR_ORANGE);
	print_event("Meio-dia solar", st.at[SOLAR_NOON], "Sol no ponto mais alto", COLOR_ORANGE);
	print_event("Hora dourada começa", st.at[GOLDEN_HOUR], "luz suave do fim da tarde", COLOR_ORANGE);
	print_event("Pôr do sol começa", st.at[SUNSET_START], "o disco toca o horizonte", COLOR_ORANGE);
	print_event("Pôr do sol", st.at[SUNSET], "o Sol desaparece", COLOR_ORANGE);
	print_event("Fim do crepúsculo civil", st.at[DUSK], "já é preciso iluminação externa", COLOR_ORANGE);
	print_event("Fim do crepúsculo náutico", st.at[NAUTICAL_DUSK], "o horizonte deixa de aparecer", COLOR_ORANGE);
	print_event("Noite astronômica", st.at[NIGHT], "céu completamente escuro", COLOR_ORANGE);

	printf("\n%s%sLUA%s\n", color(COLOR_BOLD), color(COLOR_OLIVE), color(COLOR_RESET));
	print_event("Nascer da lua", mt.rise, "a Lua aparece no horizonte", COLOR_OLIVE);
	print_event("Pôr da lua", mt.set, "a Lua desaparece no horizonte", COLOR_OLIVE);
	if (mt.always_up) printf("  A Lua permanece acima do horizonte neste dia.\n");
	if (mt.always_down) printf("  A Lua permanece abaixo do horizonte neste dia.\n");
	return 0;
}
