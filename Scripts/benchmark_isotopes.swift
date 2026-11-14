import CDKSwiftNativePort
import Foundation

// Standalone public-API probe. Build/link against the package library, then run
// in a fresh process so the first measurement includes loading the data bundle.
let coldStart = Date()
let carbon13 = try CDKIsotopes.mass(element: "C", massNumber: 13)
let coldSeconds = Date().timeIntervalSince(coldStart)
let iterations = 200_000
var checksum = 0.0
let warmStart = Date()
for index in 0..<iterations {
    let pair = [(6, 13), (43, 99), (92, 235)][index % 3]
    guard let isotope = try CDKIsotopes.isotope(atomicNumber: pair.0, massNumber: pair.1) else {
        throw CDKIsotopeError.dataUnavailable
    }
    checksum += isotope.exactMass
}
let output: [String: Any] = [
    "dataset": CDKIsotopes.datasetVersion, "carbon13": carbon13, "coldSeconds": coldSeconds,
    "lookups": iterations, "warmSeconds": Date().timeIntervalSince(warmStart), "checksum": checksum,
]
let data = try JSONSerialization.data(withJSONObject: output, options: [.sortedKeys, .prettyPrinted])
print(String(decoding: data, as: UTF8.self))
