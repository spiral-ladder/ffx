//! Provides extended finite field operations.
//!
//! Unlike `std.crypto.ff`, `Ffx` is not promised to be allocation-free nor constant time.

pub fn Ffx(comptime bits: comptime_int, comptime T: type) type {
    return struct {
        const Self = @This();

        const Fe = ff.Modulus(bits).Fe;

        m: ff.Modulus(bits),

        pub fn init(m: ff.Modulus(bits)) Self {
            return .{ .m = m };
        }

        // Computes the multiplicative inverse of a field element using Fermat's Little Theorem.
        // For a field element a, a^(-1) = a^(p-2) mod p where p is the field modulus.
        pub fn inv(self: Self, a: Fe) !Fe {
            const p_minus_2: usize = try self.m.v.toPrimitive(T) - 2;

            // Use binary exponentiation to compute a^(p-2)
            var result = self.m.one();
            var base = a;
            var exp = p_minus_2;

            while (exp > 0) {
                if (exp & 1 == 1) {
                    result = self.m.mul(result, base);
                }
                base = self.m.mul(base, base);
                exp >>= 1;
            }

            return result;
        }
    };
}

test "field element inverse" {
    const bits = 8;
    const M = ff.Modulus(bits);
    const T = u8;
    const m = try M.fromPrimitive(T, 17);
    const f = Ffx(bits, T).init(m);

    // Test inverses of all non-zero elements
    for (1..17) |i| {
        const a = M.Fe.fromPrimitive(T, m, @intCast(i)) catch unreachable;
        const a_inv = try f.inv(a);
        const prod = m.mul(a, a_inv);
        try std.testing.expectEqual(m.one(), prod);
    }

    // Test specific known inverses
    const test_cases = [_]struct { actual: u8, expected: u8 }{
        .{ .actual = 2, .expected = 9 }, // 2 * 9 = 18 ≡ 1 mod 17
        .{ .actual = 3, .expected = 6 }, // 3 * 6 = 18 ≡ 1 mod 17
        .{ .actual = 4, .expected = 13 }, // 4 * 13 = 52 ≡ 1 mod 17
        .{ .actual = 5, .expected = 7 }, // 5 * 7 = 35 ≡ 1 mod 17
        .{ .actual = 6, .expected = 3 }, // 6 * 3 = 18 ≡ 1 mod 17
        .{ .actual = 7, .expected = 5 }, // 7 * 5 = 35 ≡ 1 mod 17
        .{ .actual = 8, .expected = 15 }, // 8 * 15 = 120 ≡ 1 mod 17
        .{ .actual = 9, .expected = 2 }, // 9 * 2 = 18 ≡ 1 mod 17
    };

    for (test_cases) |case| {
        const a = M.Fe.fromPrimitive(T, m, case.actual) catch unreachable;
        const a_inv = try f.inv(a);
        const expected = M.Fe.fromPrimitive(T, m, case.expected) catch unreachable;
        try std.testing.expectEqual(expected, a_inv);
    }
}

const std = @import("std");

const ArrayList = std.ArrayList;
const Allocator = std.mem.Allocator;
const testing = std.testing;

const ff = std.crypto.ff;
