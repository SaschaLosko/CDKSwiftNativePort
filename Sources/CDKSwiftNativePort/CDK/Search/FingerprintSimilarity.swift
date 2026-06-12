import Foundation

public enum CDKFingerprintSimilarity {
    public enum SimilarityError: Error, Equatable {
        case mismatchedBitFingerprintSize
        case mismatchedVectorSize
        case emptyFingerprints
    }

    public static func tanimoto(
        _ first: CDKBitFingerprint,
        _ second: CDKBitFingerprint
    ) throws -> Double {
        guard first.size == second.size else {
            throw SimilarityError.mismatchedBitFingerprintSize
        }
        let intersection = first.bits.intersection(second.bits).count
        let union = first.bits.count + second.bits.count - intersection
        guard union > 0 else {
            throw SimilarityError.emptyFingerprints
        }
        return Double(intersection) / Double(union)
    }

    public static func tanimoto(
        _ first: CDKCountFingerprint,
        _ second: CDKCountFingerprint
    ) throws -> Double {
        try countTanimotoMethod2(first, second)
    }

    public static func tanimoto(
        _ first: [Double],
        _ second: [Double]
    ) throws -> Double {
        guard first.count == second.count else {
            throw SimilarityError.mismatchedVectorSize
        }

        var dot = 0.0
        var firstSquared = 0.0
        var secondSquared = 0.0
        for index in first.indices {
            dot += first[index] * second[index]
            firstSquared += first[index] * first[index]
            secondSquared += second[index] * second[index]
        }

        let union = firstSquared + secondSquared - dot
        guard union != 0 else {
            throw SimilarityError.emptyFingerprints
        }
        return dot / union
    }

    /// Count-fingerprint Tanimoto using the Steffen et al. dot-product method
    /// exposed as `Tanimoto.method1` in Java CDK.
    public static func countTanimotoMethod1(
        _ first: CDKCountFingerprint,
        _ second: CDKCountFingerprint
    ) throws -> Double {
        var dot: Int64 = 0
        var firstSquared: Int64 = 0
        var secondSquared: Int64 = 0

        for (hash, firstCount) in first.countsByHash {
            if let secondCount = second.countsByHash[hash] {
                dot += Int64(firstCount) * Int64(secondCount)
            }
            firstSquared += Int64(firstCount) * Int64(firstCount)
        }
        for secondCount in second.countsByHash.values {
            secondSquared += Int64(secondCount) * Int64(secondCount)
        }

        let union = firstSquared + secondSquared - dot
        guard union != 0 else {
            throw SimilarityError.emptyFingerprints
        }
        return Double(dot) / Double(union)
    }

    /// Count-fingerprint Tanimoto using the Grant et al. min/max method exposed
    /// as `Tanimoto.method2` in Java CDK.
    public static func countTanimotoMethod2(
        _ first: CDKCountFingerprint,
        _ second: CDKCountFingerprint
    ) throws -> Double {
        var minSum = 0
        var maxSum = 0
        let hashes = Set(first.countsByHash.keys).union(second.countsByHash.keys)
        for hash in hashes {
            let firstCount = first.countsByHash[hash] ?? 0
            let secondCount = second.countsByHash[hash] ?? 0
            minSum += min(firstCount, secondCount)
            maxSum += max(firstCount, secondCount)
        }
        guard maxSum != 0 else {
            throw SimilarityError.emptyFingerprints
        }
        return Double(minSum) / Double(maxSum)
    }

    public static func dice(
        _ first: CDKBitFingerprint,
        _ second: CDKBitFingerprint
    ) throws -> Double {
        guard first.size == second.size else {
            throw SimilarityError.mismatchedBitFingerprintSize
        }
        let denominator = first.bits.count + second.bits.count
        guard denominator > 0 else {
            throw SimilarityError.emptyFingerprints
        }
        let intersection = first.bits.intersection(second.bits).count
        return 2.0 * Double(intersection) / Double(denominator)
    }

    public static func cosine(
        _ first: CDKBitFingerprint,
        _ second: CDKBitFingerprint
    ) throws -> Double {
        guard first.size == second.size else {
            throw SimilarityError.mismatchedBitFingerprintSize
        }
        guard !first.bits.isEmpty, !second.bits.isEmpty else {
            throw SimilarityError.emptyFingerprints
        }
        let intersection = first.bits.intersection(second.bits).count
        return Double(intersection) / sqrt(Double(first.bits.count * second.bits.count))
    }
}
