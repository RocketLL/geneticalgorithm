package main

import (
	"fmt"
	"image"
	"image/color"
	"image/draw"
)

func DrawCircle(img draw.Image, x0, y0, r int, c color.Color) {
	for y := -r; y <= r; y++ {
		for x := -r; x <= r; x++ {
			if x*x+y*y <= r*r {
				img.Set(x0+x, y0+y, c)
			}
		}
	}
}

func UnTransparent(img draw.Image) {
	white := color.RGBA{255, 255, 255, 255}
	draw.Draw(img, img.Bounds(), &image.Uniform{white}, image.ZP, draw.Src)
}

//func DrawCircle(img draw.Image, x0, y0, r int, c color.Color) {
//	x, y, dx, dy := r-1, 0, 1, 1
//	err := dx - (r * 2)
//
//	for x > y {
//		img.Set(x0+x, y0+y, c)
//		img.Set(x0+y, y0+x, c)
//		img.Set(x0-y, y0+x, c)
//		img.Set(x0-x, y0+y, c)
//		img.Set(x0-x, y0-y, c)
//		img.Set(x0-y, y0-x, c)
//		img.Set(x0+y, y0-x, c)
//		img.Set(x0+x, y0-y, c)
//
//		if err <= 0 {
//			y++
//			err += dy
//			dy += 2
//		}
//		if err > 0 {
//			x--
//			dx += 2
//			err += dx - (r * 2)
//		}
//	}
//}

//func CompareImage(img1, img2 *image.RGBA) (int64, error) {
//	if img1.Bounds() != img2.Bounds() {
//		return 0, fmt.Errorf("image bounds not equal: %+v, %+v", img1.Bounds(), img2.Bounds())
//	}
//
//	accumError := int64(0)
//
//	for i := 0; i < len(img1.Pix); i++ {
//		accumError += int64(sqDiffUInt8(img1.Pix[i], img2.Pix[i]))
//	}
//
//	return int64(math.Sqrt(float64(accumError))), nil
//}

// func CompareImage(img1, img2 *image.RGBA) (int64, error) {
// 	if img1.Bounds() != img2.Bounds() {
// 		return 0, fmt.Errorf("image bounds not equal: %+v, %+v", img1.Bounds(), img2.Bounds())
// 	}

// 	accumError := int64(0)

// 	for i := 0; i < len(img1.Pix); i++ {
// 		accumError += int64(sqDiffUInt8(img1.Pix[i], img2.Pix[i]))
// 	}

// 	return int64(math.Sqrt(float64(accumError))), nil
// }

func CompareImage(img1, img2 *image.RGBA) (float64, error) {
	b := img1.Bounds()
	if b != img2.Bounds() {
		return 0, fmt.Errorf("image bounds not equal: %+v, %+v", img1.Bounds(), img2.Bounds())
	}

	accumError := int64(0)

	for i := 0; i < len(img1.Pix); i++ {
		if i%4 != 3 {
			err := int64(int(img1.Pix[i]) - int(img2.Pix[i]))
			accumError += err * err
		}
	}

	mse := float64(accumError) / float64(3*b.Max.X*b.Max.Y)
	return mse, nil
}
