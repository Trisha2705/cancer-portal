<?php
session_start();

require_once __DIR__ . '/../includes/db.php';

// Ensure the user is logged in.
if (!isset($_SESSION['user_id']) || $_SESSION['role'] !== 'professional') {
    header('Location: ../login.php');
    exit;
}

$user_id = (int) $_SESSION['user_id'];

$prediction_id = filter_input(
    INPUT_GET,
    'prediction_id',
    FILTER_VALIDATE_INT
);

if (!$prediction_id || $prediction_id < 1) {
    die('Invalid prediction ID.');
}

$reviewMessage = '';
$reviewError = '';

function h($value) {
    return htmlspecialchars((string) $value, ENT_QUOTES, 'UTF-8');
}

try {

    // Verify that this professional is an oncologist.
    $stmt = $pdo->prepare("
        SELECT professional_id, first_name, last_name
        FROM professionals
        WHERE user_id = ?
          AND professional_type = 'oncologist'
    ");

    $stmt->execute([$user_id]);
    $professional = $stmt->fetch(PDO::FETCH_ASSOC);

    if (!$professional) {
        die('Oncologist access required.');
    }

    // Retrieve prediction and verify patient assignment.
    $stmt = $pdo->prepare("
        SELECT
            p.prediction_id,
            p.patient_id,
            p.model_used,
            p.created_at,
            pat.first_name AS patient_first_name,
            pat.last_name AS patient_last_name
        FROM predictions p
        INNER JOIN patients pat
            ON pat.patient_id = p.patient_id
        INNER JOIN oncologist_patients op
            ON op.patient_id = p.patient_id
        WHERE p.prediction_id = ?
          AND op.user_id = ?
          AND op.assignment_status = 'active'
    ");

    $stmt->execute([$prediction_id, $user_id]);
    $prediction = $stmt->fetch(PDO::FETCH_ASSOC);

    if (!$prediction) {
        die('Prediction not found or access denied.');
    }

    $patient_id = (int) $prediction['patient_id'];

    // Retrieve prediction results.
    $stmt = $pdo->prepare("
        SELECT
            ct.cancer_name,
            ct.cancer_code,
            pr.prediction_probability
        FROM prediction_results pr
        INNER JOIN cancer_types ct
            ON ct.cancer_type_id = pr.cancer_type_id
        WHERE pr.prediction_id = ?
        ORDER BY pr.prediction_probability DESC
    ");

    $stmt->execute([$prediction_id]);
    $predictionResults = $stmt->fetchAll(PDO::FETCH_ASSOC);

    // Handle review submission.
    if ($_SERVER['REQUEST_METHOD'] === 'POST'
        && isset($_POST['submit_review'])) {

        $agreement_status = $_POST['agreement_status'] ?? '';
        $review_reasoning = trim($_POST['review_reasoning'] ?? '');

        $validStatuses = [
            'agree',
            'partially_agree',
            'disagree'
        ];

        if (!in_array($agreement_status, $validStatuses, true)) {
            $reviewError = 'Please select an agreement status.';
        } elseif ($review_reasoning === '') {
            $reviewError = 'Please enter your review reasoning.';
        } else {

            $stmt = $pdo->prepare("
                INSERT INTO oncologist_reviews (
                    user_id,
                    patient_id,
                    prediction_id,
                    review_title,
                    agreement_status,
                    review_notes,
                    review_reasoning,
                    review_status,
                    reviewed_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, 'completed', NOW())
            ");

            $title = 'AI Prediction Review #' . $prediction_id;

            $stmt->execute([
                $user_id,
                $patient_id,
                $prediction_id,
                $title,
                $agreement_status,
                $review_reasoning,
                $review_reasoning
            ]);

            header(
                'Location: prediction-review.php?prediction_id='
                . $prediction_id . '&saved=1'
            );
            exit;
        }
    }

    if (isset($_GET['saved']) && $_GET['saved'] === '1') {
        $reviewMessage = 'AI prediction review saved successfully.';
    }

    // Retrieve previous reviews for this prediction.
    $stmt = $pdo->prepare("
        SELECT
            agreement_status,
            review_reasoning,
            created_at
        FROM oncologist_reviews
        WHERE prediction_id = ?
          AND user_id = ?
        ORDER BY created_at DESC
    ");

    $stmt->execute([$prediction_id, $user_id]);
    $previousReviews = $stmt->fetchAll(PDO::FETCH_ASSOC);

} catch (PDOException $e) {

    error_log('Prediction review page error: ' . $e->getMessage());
    die('Unable to load the prediction review page.');
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AI Prediction Review | CancerInsight</title>
    <link rel="stylesheet" href="../style/app.css">
</head>

<body>

<main class="page-container">

    <h1>AI Prediction Review</h1>

    <p>
        Patient:
        <strong>
            <?= h($prediction['patient_first_name'] . ' ' . $prediction['patient_last_name']) ?>
        </strong>
    </p>

    <p>
        <strong>Prediction ID:</strong>
        <?= (int) $prediction['prediction_id'] ?>
    </p>

    <p>
        <strong>Model:</strong>
        <?= h($prediction['model_used']) ?>
    </p>

    <p>
        <strong>Generated:</strong>
        <?= h($prediction['created_at']) ?>
    </p>

    <p>
        <a href="prediction-review.php?prediction_id=<?= (int) $prediction['prediction_id'] ?>">
            Review this AI prediction
        </a>
    </p>

    <section class="details-section">

        <h2>AI Prediction Results</h2>

        <?php if (empty($predictionResults)): ?>

            <p>No prediction results are available.</p>

        <?php else: ?>

            <?php foreach ($predictionResults as $result): ?>

                <article class="record-card">

                    <h3>
                        <?= h($result['cancer_name']) ?>
                        (<?= h($result['cancer_code']) ?>)
                    </h3>

                    <p>
                        <strong>Probability:</strong>
                        <?= number_format(
                            (float) $result['prediction_probability'] * 100,
                            2
                        ) ?>%
                    </p>

                </article>

            <?php endforeach; ?>

        <?php endif; ?>

    </section>

    <section class="details-section">

        <h2>Clinician Review</h2>

        <p>
            Record your professional review of the AI-generated output.
        </p>

        <?php if ($reviewMessage): ?>
            <p class="success-message">
                <?= h($reviewMessage) ?>
            </p>
        <?php endif; ?>

        <?php if ($reviewError): ?>
            <p class="error-message">
                <?= h($reviewError) ?>
            </p>
        <?php endif; ?>

        <form method="POST">

            <div>

                <label for="agreement_status">
                    Agreement with AI prediction
                </label>

                <select
                    name="agreement_status"
                    id="agreement_status"
                    required
                >
                    <option value="">Select an option</option>
                    <option value="agree">Agree</option>
                    <option value="partially_agree">
                        Partially agree
                    </option>
                    <option value="disagree">Disagree</option>
                </select>

            </div>

            <div>

                <label for="review_reasoning">
                    Review reasoning
                </label>

                <textarea
                    name="review_reasoning"
                    id="review_reasoning"
                    rows="6"
                    required
                    placeholder="Explain your review of the prediction..."
                ></textarea>

            </div>

            <button
                type="submit"
                name="submit_review"
                value="1"
                class="ci-button"
            >
                Save Review
            </button>

        </form>

    </section>

    <section class="details-section">

        <h2>Previous Reviews</h2>

        <?php if (empty($previousReviews)): ?>

            <p>No reviews have been recorded for this prediction.</p>

        <?php else: ?>

            <?php foreach ($previousReviews as $review): ?>

                <article class="record-card">

                    <p>
                        <strong>Agreement:</strong>
                        <?= h(str_replace(
                            '_',
                            ' ',
                            ucfirst($review['agreement_status'])
                        )) ?>
                    </p>

                    <p>
                        <strong>Reasoning:</strong>
                    </p>

                    <p>
                        <?= nl2br(h($review['review_reasoning'])) ?>
                    </p>

                    <p>
                        <strong>Date:</strong>
                        <?= h($review['created_at']) ?>
                    </p>

                </article>

            <?php endforeach; ?>

        <?php endif; ?>

    </section>

    <p>
        <a href="patients.php">Back to Patients</a>
    </p>

</main>

</body>
</html>
